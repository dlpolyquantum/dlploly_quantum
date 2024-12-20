#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --time 48:00:00
#SBATCH --ntasks-per-node 1
#SBATCH --mail-user=nlondon@umkc.edu
#SBATCH --mail-type=FAIL


module purge
module load intel-mpi
module load fftw
module list

set -x 
echo JOB STARTED AT:
 date

[ -e Average_Intra.d ] && rm Average_Intra.d 
[ -e Average_Angle.d ] && rm Average_Angle.d 
[ -e Average_RDF.d ] && rm Average_RDF.d 

export EXE=/jet/home/nl478/programs/rdfAvg/rdfAvg
#mpiexec.hydra -bootstrap sge $EXE
$EXE

echo JOB FINISHED AT:
 date
echo ALL DONE!

