#!/bin/bash
mkdir ref
mkdir ref/run
cp base_scripts/run_dlpoly_fqcmd.sh ref/run/
cp base_scripts/run_avg.sh ref/run/
cd ref/run
mkdir base_input
cp ../../FQCMD base_input/
cp ../../refCONTROL base_input/CONTROL
cp ../../avgCONTROL CONTROL
cp ../../zero-pot/*.table base_input/ 
cp ../../zero-pot/FIELD base_input/FIELD
for k in $(seq -f "%03g" 2 50)
do
	mkdir traj${k}
	cp run_dlpoly_fqcmd.sh traj${k}/
	cp ../../config/CONFIG${k} traj${k}/CONFIG 
	cd traj${k}
	sbatch -J ref-${k} run_dlpoly_bridges2.sh
	cd ../
done	
cd ../../	
