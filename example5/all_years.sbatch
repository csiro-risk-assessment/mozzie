#!/bin/bash -l
#SBATCH --job-name=calc_qm
##### NOTE: need to check the following time is appropriate
#SBATCH --time=3:00:00
#SBATCH --nodes=1
##### NOTE: assuming 10 years, run in parallel
#SBATCH --ntasks=10
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.wilkins@csiro.au
##SBATCH --export=NONE

#### Note: to check the CPU usage of your job (should be 100%), do something like
#### ssh -t c201 top -u wil04q -d1 -c -n 10

module load SC slurm python/3.9.4 parallel

##### NOTE: assuming 10 years, run in parallel
years="2010 2011 2012 2013 2014 2015 2016 2017 2018 2019"
parallel -j10 --xapply 'python yearly_qm.py ${1} &> {1}.log' ::: ${years}
