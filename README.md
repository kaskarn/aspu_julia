# aspu_julia
aspu implemented in julia
```

Adaptive Sum of Powered Scores (aSPU)
        Reference: doi:10.1002/gepi.21931
        Contact: Antoine Baldassari baldassa@email.unc.edu

        Computes aggregate p-values from multiple traits, using GWAS results.

        Input is a comma-separated file with a header, with SNP names in the
        first column, and each subsequent column containing summary z-scores.

        Output is comma-separated, with columns for SNP name, aspu p-value,
        and best-powered gamma value

        The SLURM log file will be named aspu_julia_<timestamp>.out


        USAGE

        bash aspu_slurm.sh --ncpu int --filein filename --logB int ...

        Example: bash aspu_slurm.sh --ncpu 10 --filein myfile.csv --logB 6
                 --fileout yaymath.txt --slurm "--begin=now+1hour"
                 --name myaspuresults


        ARGUMENTS

        --name dirname [default: aspu_run_date_time]
        Name of directory where results will be stored.

        --filein filename (REQUIRED)
        path to input file

        --logB integer between 1 and 9 (REQUIRED)
        log10 of the number of Monte-Carlo simulations
        For example, --logB 8 will allow the algorithm to run up to 10^8
        simulations, allowing p-values as small as 10^-8

        --ncpu integer (REQUIRED)
        The number of CPUs to request from the cluster
        (more than 20 is likely overkill)

        --fileout filename
        Name of output file.
        The default is aspu_results_<timestamp>_1E<logB>.csv

        --incov filename
        Specify own (comma-separated) file containing covariance matrix
        of z-scores. The default is to calculate if from input results

        --outcov filename
        Specify name of output z-scores covariance matrix

        --norun
        Stop program after setup phase. Used to check for mistakes and
        verify arguments before commiting to full run (check log).

        --slurm string
        SLURM options to pass to sbatch. Must be enclosed in " or ' quotes
        and cannot conflict with options -n, --mem-per-cpu, --cpus-per-task,
        and -o

        --pows
        Specify powers to be used. Integer powers must be ranges, or
        comma-separated integers.

        The integer 9 is a special value reserved for gamma=infty (MinP)
        Example: 1:5,9
        Default is 1:9

        --testnow
        Run aspu in the local session instead of submitting using SLURM
        For testing.

```
