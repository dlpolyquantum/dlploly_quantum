#!/bin/bash

#Bash script for automatically setting up and running 
#IBI iterations with multiple trajectories including
#depencies to wait for each iteration to finish before
#starting next one
#
#Author: Dil Limbu, 2024

int_itr=1
max_itr=30
ntraj=50

for (( i=int_itr; i<=max_itr; i++ )); 
do
let j=i-1

echo "Starting Iteration $i $j ..."

if [ -d iter${i} ]; then
  echo iter${i} exists, change iteration number at first!!
  exit 1
elif [ ! -d iter${j} ]; then
  echo iter${j} does not exist, change iteration number at first!!
  exit 1
else
  echo "running simulation for iter${i}"
  echo "copying necessary files from iter${j}"
fi
mkdir iter${i}
mkdir iter${i}/run
cp base_scripts/run_avg.sh iter${i}/run/
cp base_scripts/run_copyfiles.sh iter${i}/run/
cp base_scripts/run_IBI.sh iter${i}/
cd iter${i}/run
mkdir base_input
cp ../../FQCMD base_input/
cp ../../iterCONTROL base_input/CONTROL
cp ../../avgCONTROL CONTROL

for k in $(seq -f "%03g" 1 ${ntraj})
do
	mkdir traj${k}
	cp ../../base_scripts/run_dlpoly_fqcmd.sh traj${k}/
	cp ../../config/CONFIG${k} traj${k}/CONFIG 
done

if [ "$i" -eq "$int_itr" ]; then
 run_copy_id=$( sbatch -J copy${i} --parsable run_copyfiles.sh $j )
 echo "Submitted Copy Files for iteration $i with Job ID $run_copy_id"
else
 run_copy_id=$( sbatch -J copy${i} --parsable --dependency=afterok:${ibi_pre_job} run_copyfiles.sh $j )
 echo "Submitted Copy Files for iteration $i with Job ID $run_copy_id after $ibi_pre_job"
fi
echo 'run_copy_id' $run_copy_id

job_ids=()
for k in $(seq -f "%03g" 1 ${ntraj})
do
   cd traj${k}
   job_id=$(sbatch -J it${i}-${k} --parsable --dependency=afterok:${run_copy_id} run_dlpoly_fqcmd.sh)
   job_ids+=($job_id)
   echo $job_id
   cd ../
done

# Set dependency for post-analysis job to run after all 20 jobs are complete
dependency_list=$(IFS=: ; echo "${job_ids[*]}")
# Submit the analysis script with a dependency on the completion of the setup job
ave_job_id=$(sbatch -J avgRDF${i} --parsable --dependency=afterok:${dependency_list} run_avg.sh)
echo "Submitted avgRDF for iteration $i with Job ID $ave_job_id after $dependency_list"
cd ..

ibi_job_id=$(sbatch -J runIBI${i} --parsable --dependency=afterok:${ave_job_id} run_IBI.sh)
echo "Submitted IBI run for iteration $i with Job ID $ibi_job_id" after $ave_job_id
ibi_pre_job=$ibi_job_id
cd ..

done
