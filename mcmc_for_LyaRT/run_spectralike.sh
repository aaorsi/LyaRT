#!/bin/tcsh
#$ -pe orte-20 60
#$ -cwd
#$ -S /bin/tcsh

set mcmc_dir = '/home/CEFCA/aaorsi/mcmc_for_LyaRT'

mpirun -np 60 python $mcmc_dir/spectralike.py

