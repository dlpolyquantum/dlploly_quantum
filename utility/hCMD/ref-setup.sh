#!/bin/bash
mkdir ref
mkdir ref/run
cp base_scripts/run_dlpoly_fqcmd.sh ref/run/
cp base_scripts/run_avg.sh ref/run/
cd ref/run
mkdir base_input
cp ../../FQCMD base_input/
cp ../../refCONTROL base_input/CONTROL
cp ../../refFIELD base_input/FIELD
cp ../../avgCONTROL CONTROL
cp ../../zero-pot/*.table base_input/ 
for k in $(seq -f "%03g" 1 50)
do
	mkdir traj${k}
	cp run_dlpoly_fqcmd.sh traj${k}/
	cp ../../config/CONFIG${k} traj${k}/CONFIG 
	cd traj${k}
	sbatch -J ref-${k} run_dlpoly_fqcmd.sh
	cd ../
done	
cd ../../	
