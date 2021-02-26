#!/bin/bash
#QUEUE=zeus_short_q
QUEUE=zeus_all_q
#QUEUE=zeus_long_q
for p in $(seq 1.25 1 1.25)
#for p in $(seq 0.75 0.25 1.25)
do
	for (( f=1; f<=4; f++ ))
	do
		for (( c=2; c<=2; c++ ))
		do
			echo "________________________________________________"
			echo "Submitting to Zeus: fuelIdx=$f, iter=$c, psf=$p."
			qsub -v fid=$f,iter=$c,psf=$p -N I${c}_F${f}_P${p} -q $QUEUE run_sim_z_1p.sh
			sleep 0.2
		done	
	done
done