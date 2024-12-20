#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --time 48:00:00
#SBATCH --ntasks-per-node 1
#SBATCH --mail-user=nlondon@umkc.edu
#SBATCH --mail-type=ALL

module purge
module load intel-mpi
module load fftw
module load anaconda3
module list

echo JOB STARTED AT:
 date

#export EXE=/jet/home/nl478/programs/rdfAvg/rdfAvg
#mpiexec.hydra -bootstrap sge $EXE
#$EXE
python ../../../../IBI_water.py

echo JOB FINISHED AT:
 date
echo ALL DONE!

