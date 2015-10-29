#!/bin/bash -l
#SBATCH -D /home/jolyang/Documents/Github/imputeR
#SBATCH -J teodata
#SBATCH -o /home/jolyang/Documents/Github/imputeR/slurm-log/out-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/imputeR/slurm-log/error-%j.txt

R --no-save < profiling/cjdata/1.A.1_loading_recoding.R
#sbatch -p bigmemh --mem 160000 --ntasks=20 slurm-script/run-teodata.sh