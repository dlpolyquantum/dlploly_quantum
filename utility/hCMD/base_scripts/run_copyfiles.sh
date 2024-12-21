#!/bin/bash -l
#SBATCH --job-name=copyFiles
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH --time 00:10:00
#SBATCH --ntasks-per-node 1


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

sleep 5

echo JOB FINISHED AT:
 date
echo ALL DONE!

