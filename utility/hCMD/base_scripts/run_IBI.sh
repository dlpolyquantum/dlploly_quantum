#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --time 00:10:00
#SBATCH --ntasks-per-node 1

module purge
module load intel-mpi
module load fftw
module load anaconda3
module list

echo JOB STARTED AT:
 date

python ../IBI_litfsi.py

echo JOB FINISHED AT:
 date
echo ALL DONE!

