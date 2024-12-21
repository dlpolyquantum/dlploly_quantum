#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --time 48:00:00
#SBATCH --ntasks-per-node 32
#SBATCH --mem=20000mb

module purge
module load intel-mpi
module list

echo JOB STARTED AT:
 date

#Set to filepath of DLPOLY.X
#export EXE=
export EXE=/jet/home/limbu/Softwares/dlpoly_quantum_v2/execute/DLPOLY.X

mpiexec.hydra -bootstrap sge $EXE

echo JOB FINISHED AT:
 date
echo ALL DONE!

