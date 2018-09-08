#!/bin/bash
homedir=`pwd`
script=$(readlink -f "$0")
scriptpath=$(dirname "$script")

if [ $# -eq 0 ]; then
	more $scriptpath/README.md
  exit
fi

#Parse input parameters
while [[ $# -gt 0 ]]; do
  if [ "$1" == "--filein" ]; then
    filein=`readlink -f $2`
    tojulia="$tojulia $1 $filein "
    shift 2
    continue
  fi 
  if [ "$1" == "--incov" ]; then
    incov=`readlink -f $2`
    tojulia="$tojulia $1 $incov "
    shift 2
    continue
  fi 
	tojulia="$tojulia $1"
	[[ "$1" == "--ncpu" ]] && ncpu=$2
	[[ "$1" == "--norun" ]] && norun="_norun"
	[[ "$1" == "--logB" ]] && logB=$2
	[[ "$1" == "--slurm" ]] && slurmopts=$2
  [[ "$1" == "--testnow" ]] && testnow="true"
  [[ "$1" == "--name" ]] && name=$2
	shift 
done

#Detect errors
[[ -z $filein ]] && err="ERROR: Option --filein is required."
[[ -z $logB ]] && err="ERROR: Option --logB is required.\n$err"
[[ -z $ncpu ]] && err="ERROR: Option --ncpu is required.\n$err"
echo -e "$err"
[[ ! -z $err ]] && exit

#Create results directory
[[ -z $name ]] && name="aspu_run_`date +%Y_%m_%d_%Hh_%Mm_%Ss`"
[[ -d $name ]] || mkdir $name
cd $name

#Make covariance matrix
if [ -z "$incov" ]; then
  echo "Creating covariance matrix from $filein..."
  tmr=`mktemp`
  covnam="covmat_`basename $filein`"
  echo "
library(data.table)
fin<-fread(\"$filein\")
keep<-fin[,do.call(\"pmax\",lapply(.SD[,-1], abs))] < (-qnorm(0.5e-5,lower.tail = T))
rmat<-cor(fin[keep,-1])
write.table(rmat,\"$covnam\",col.names=F,row.names=F,sep=\",\")
" > $tmr
  Rscript $tmr
  rm $tmr
  tojulia="$tojulia --incov `readlink -f $covnam`"
fi

#Set exe paths
jexec="julia"
juliacall="$jexec $scriptpath/julia/aspu_io.jl"

#Install packages as needed
[[ -d ~/.julia/packages ]] || mkdir -p ~/.julia/packages
addpkg=`comm -13 <(ls ~/.julia/packages/) <(echo -e "ClusterManagers\nCSV\nDistributions")`
[[ -z $addpkg ]] || echo "Installing missing Julia Packages. This will take a few minutes.\nPackage(s): $addpkg ..."
[[ -z $addpkg ]] || $jexec -e "using Pkg; [Pkg.add(i) for i = [`sed 's/^\|$/"/g' <(echo $addpkg) | sed 's/ /" "/g'`]]"

acct=`sacctmgr list User format="DefaultAccount%30,User" | grep $USER | xargs | cut -f1 -d' '`
mem="5GB"
# save command
echo sbatch -o "aspu_julia${norun}_`date +%Y_%m_%d_%Hh_%Mm_%Ss`.out" -A $acct -n $ncpu --cpus-per-task 1 -N 2-$ncpu --time=7-0 --mem-per-cpu=$mem $slurmopts --wrap="$juliacall $tojulia" > aspu_cmd.txt

[[ -z $testnow ]] && sbatch -o "aspu_julia${norun}_`date +%Y_%m_%d_%Hh_%Mm_%Ss`.out" -A $acct -n $ncpu --cpus-per-task 1 -N 2-$ncpu --time=7-0 --mem-per-cpu=$mem $slurmopts --wrap="$juliacall $tojulia"
[[ -z $testnow ]] || $juliacall $tojulia


cd $homedir


