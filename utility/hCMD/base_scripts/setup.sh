#!/bin/bash
for i in $(seq -f "%03g" 2 50)
do
	echo $i
	mkdir traj${i}
	cp base_input/* traj${i}/
	cp run_dlpoly_bridges2.sh traj${i}/
#	mv config/CONFIG${i} traj${i}/CONFIG
	cp ../../ref/run/traj${i}/CONFIG traj${i}/CONFIG
done

