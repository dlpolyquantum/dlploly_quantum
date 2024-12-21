#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --time 1:00:00
#SBATCH --ntasks-per-node 1

module purge
module load intel-mpi
module load fftw
module list

echo JOB STARTED AT:
 date

export EXE=/jet/home/nl478/programs/correlation/correlation
#mpiexec.hydra -bootstrap sge $EXE
$EXE

echo JOB FINISHED AT:
 date
echo ALL DONE!

