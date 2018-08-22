#!/bin/bash
homedir=`pwd`
if [ $# -eq 0 ]; then
	echo -e "\n\n\tAdaptive Sum of Powered Scores (aSPU)\n\tReference: doi:10.1002/gepi.21931\n\tContact: Antoine Baldassari baldassa@email.unc.edu"
	echo -e "\n\tComputes aggregate p-values from multiple traits, using GWAS results."
	echo -e "\n\tInput is a comma-separated file with a header, with SNP names in the"
	echo -e "\tfirst column, and each subsequent column containing summary z-scores.\n"
	echo -e "\tOutput is comma-separated, with columns for SNP name, aspu p-value,"
	echo -e "\tand best-powered gamma value\n"
	echo -e "\tThe SLURM log file will be named aspu_julia_<timestamp>.out\n"
	echo -e "\n\tUSAGE\n"
	echo -e "\tbash aspu_slurm.sh --ncpu int --filein filename --logB int ...\n"
	echo -e "\tExample: bash aspu_slurm.sh --ncpu 10 --filein myfile.csv --logB 6"
	echo -e "\t         --fileout yaymath.txt --slurm \"--begin=now+1hour\""
	echo -e "\t         --name myaspuresults\n"
	echo -e "\n\tARGUMENTS\n"
  echo -e "\t--name dirname [default: aspu_run_date_time]\n\tName of directory where results will be stored.\n"
	echo -e "\t--filein filename (REQUIRED)\n\tpath to input file\n"
	echo -e "\t--logB integer between 1 and 9 (REQUIRED)\n\tlog10 of the number of Monte-Carlo simulations"
	echo -e "\tFor example, --logB 8 will allow the algorithm to run up to 10^8"
	echo -e "\tsimulations, allowing p-values as small as 10^-8\n"
	echo -e "\t--ncpu integer (REQUIRED)\n\tThe number of CPUs to request from the cluster"
	echo -e "\t(more than 20 is likely overkill)\n"
	echo -e "\t--fileout filename\n\tName of output file.\n\tThe default is aspu_results_<timestamp>_1E<logB>.csv\n"
	echo -e "\t--incov filename\n\tSpecify own (comma-separated) file containing covariance matrix"
	echo -e "\tof z-scores. The default is to calculate if from input results\n"
	echo -e "\t--outcov filename\n\tSpecify name of output z-scores covariance matrix\n"
	echo -e "\t--norun\n\tStop program after setup phase. Used to check for mistakes and"
	echo -e "\tverify arguments before commiting to full run (check log).\n"
	echo -e "\t--slurm string\n\tSLURM options to pass to sbatch. Must be enclosed in \" or ' quotes"
	echo -e "\tand cannot conflict with options -n, --mem-per-cpu, --cpus-per-task,"
	echo -e "\tand -o\n"
  echo -e "\t--pows\n\tSpecify powers to be used. Integer powers must be ranges, or\n\tcomma-separated integers.\n"
  echo -e "\tThe integer 9 is a special value reserved for gamma=infty (MinP)\n\tExample: 1:5,9\n\tDefault is 1:9\n"
  echo -e "\t--testnow\n\tRun aspu in the local session instead of submitting using SLURM\n\tFor testing.\n"
	echo -e "\n"
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

[[ -z $norun ]] || ncpu=2
 
# juliacall="$jexec /proj/epi/CVDGeneNas/antoine/bin/aspu _julia/aspu_alt.jl"
juliacall="$jexec /proj/epi/CVDGeneNas/antoine/bin/aspu_julia_1.0/aspu_io.jl"

# echo $tojulia
echo sbatch -o "aspu_julia${norun}_`date +%Y_%m_%d_%Hh_%Mm_%Ss`.out" -n $ncpu --cpus-per-task 1 -N 2-$ncpu --mem-per-cpu=$mem --time=7-0 $slurmopts --wrap="$juliacall $tojulia" > aspu_log.txt
[[ -z $testnow ]] && sbatch -o "aspu_julia${norun}_`date +%Y_%m_%d_%Hh_%Mm_%Ss`.out" -n $ncpu --cpus-per-task 1 -N 2-$ncpu --mem-per-cpu=$mem --time=7-0 $slurmopts --wrap="$juliacall $tojulia"
[[ -z $testnow ]] || $juliacall $tojulia

cd $homedir


