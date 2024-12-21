#!/bin/bash -l
#SBATCH --job-name=avgRDF
#SBATCH --nodes=1
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node 1
#SBATCH --time=0-00:10:00

set -x 
echo JOB STARTED AT:
 date

[ -e Average_Intra.d ] && rm Average_Intra.d 
[ -e Average_Angle.d ] && rm Average_Angle.d 
[ -e Average_RDF.d ] && rm Average_RDF.d 

export EXE=/jet/home/limbu/Softwares/codes/programs/rdfAvg/rdfAvg
#mpiexec.hydra -bootstrap sge $EXE
$EXE

echo JOB FINISHED AT:
 date
echo ALL DONE!

