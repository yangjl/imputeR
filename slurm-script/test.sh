#!/bin/bash -l
#SBATCH -D /home/jolyang/Documents/Github/phasing_tests
#SBATCH -J phase_test
#SBATCH -o /home/jolyang/Documents/Github/phasing_tests/slurm/phase_out-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/phasing_tests/slurm/phase_error-%j.txt
#SBATCH --array=1-150

#only argument is mean (of poisson) of number of crossovers
crossovers=$1
R --no-save "--args $crossovers $SLURM_JOB_ID" < test_code.r 
