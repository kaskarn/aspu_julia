#!/bin/bash
homedir=`pwd`
script=$(readlink -f "$0")
scriptpath=$(dirname "$script")

if [ $# -eq 0 ]; then
	more $scriptpath/README.md
  exit
fi

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

[[ -z $filein ]] && err="ERROR: Option --filein is required."
[[ -z $logB ]] && err="ERROR: Option --logB is required.\n$err"
[[ -z $ncpu ]] && err="ERROR: Option --ncpu is required.\n$err"
echo -e "$err"
[[ ! -z $err ]] && exit

[[ -z $name ]] && name="aspu_run_`date +%Y_%m_%d_%Hh_%Mm_%Ss`"
[[ -d $name ]] || mkdir $name
cd $name

filein=`readlink -f $filein`
[[ -z $incov ]] || incov=`readlink -f $incov`


mem="7GB"
jexec="/proj/epi/CVDGeneNas/antoine/bin/julia-1.0.0/bin/julia"

[[ -d ~/.julia/packages ]] || mkdir -p ~/.julia/packages
# echo "Updating Julia packages... this may take a few minutes"
# $jexec -e "using Pkg; Pkg.update()"

addpkg=`comm -13 <(ls ~/.julia/packages/) <(echo -e "ClusterManagers\nCSV\nDistributions")`
[[ -z $addpkg ]] || echo "Installing missing Julia Packages. This will take a few minutes.\nPackage(s): $addpkg ..."
[[ -z $addpkg ]] || $jexec -e "using Pkg; [Pkg.add(i) for i = [`sed 's/^\|$/"/g' <(echo $addpkg) | sed 's/ /" "/g'`]]"

# [[ -z $norun ]] || ncpu=2
 
juliacall="$jexec $scriptpath/julia/aspu_io.jl"

# echo $tojulia
echo sbatch -o "aspu_julia${norun}_`date +%Y_%m_%d_%Hh_%Mm_%Ss`.out" -n $ncpu --cpus-per-task 1 -N 1-$ncpu --mem-per-cpu=$mem --time=7-0 $slurmopts --wrap="$juliacall $tojulia" > aspu_log.txt

[[ -z $testnow ]] && sbatch -o "aspu_julia${norun}_`date +%Y_%m_%d_%Hh_%Mm_%Ss`.out" -n $ncpu --cpus-per-task 1 -N 1-$ncpu --mem-per-cpu=$mem --time=7-0 $slurmopts --wrap="$juliacall $tojulia"
[[ -z $testnow ]] || $juliacall $tojulia
# [[ -z $testnow ]] || srun --pty -n $ncpu --cpus-per-task 1 -N 2-$ncpu --mem-per-cpu=$mem --time=7-0 $slurmopts $juliacall $tojulia

cd $homedir


