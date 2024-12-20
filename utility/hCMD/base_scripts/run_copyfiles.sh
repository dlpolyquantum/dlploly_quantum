#!/bin/bash -I
#SBATCH --job-name=runAVG
#SBATCH --output=myoutput%j.out
#SBATCH --error=jobError.out
#SBATCH --time=0:05:00
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=fail
#SBATCH --mail-user=nlondon@umkc.edu

echo JOB STARTED AT:
 date

j=$1

echo "Copying files from iter$j ..."

cp ../../iter${j}/*.table base_input/ 
cp ../../iter${j}/new_FIELD base_input/FIELD
cp ../../iter${j}/new_FIELD ../old_FIELD
cp ../../iter${j}/New_Potential.dat ../Previous_Potential.dat

for k in traj*
do
#	mkdir traj${k}
	cp base_input/* ${k}/
done

sleep 10

echo JOB FINISHED AT:
 date
echo ALL DONE!

