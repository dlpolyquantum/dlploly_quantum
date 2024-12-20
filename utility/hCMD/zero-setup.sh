#!/bin/bash
let i=0
echo ${j}
mkdir iter${i}
mkdir iter${i}/run
cp base_scripts/run_dlpoly_fqcmd.sh iter${i}/run/
cp base_scripts/run_avg.sh iter${i}/run/
cd iter${i}/run
mkdir base_input
cp ../../FQCMD base_input/
cp ../../iterCONTROL base_input/CONTROL
cp ../../avgCONTROL CONTROL
cp ../../zero-pot/*.table base_input/ 
cp ../../zero-pot/FIELD base_input/FIELD
cp ../../zero-pot/FIELD ../old_FIELD
cp ../../zero-pot/Potential.dat ../Previous_Potential.dat
for k in $(seq -f "%03g" 1 50)
do
	mkdir traj${k}
#	cp base_input/* traj${k}/
	cp run_dlpoly_fqcmd.sh traj${k}/
	cp ../../config/CONFIG${k} traj${k}/CONFIG 
	cd traj${k}
	sbatch -J it${i}-${k} run_dlpoly_bridges2.sh
	cd ../
done	
cd ../../	
