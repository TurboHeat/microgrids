#!/bin/bash
CPUs=18
QUEUE=zeus_all_q
#QUEUE=zeus_long_q

#for p in $(seq 0.75 0.25 1.25)
#for p in $(seq 1 0.25 1.25)
#for p in $(seq 0.75 1 0.75)
#for p in $(seq 1 1 1)
for p in $(seq 1.25 1.00 1.25)
do
	for (( f=1; f<=4; f++ ))
	do
		m=138
		M=155
		echo "______________________________________________________"
		echo "Submitting to Zeus: fuelIdx=$f, iter=[$m..$M], psf=$p."
		qsub -v fid=$f,minid=$m,maxid=$M,W=$CPUs,psf=$p -N I${m}-${M}_F${f}_P${p} -q $QUEUE -l select=1:ncpus=$CPUs -l place=scatter run_sim_z_mp.sh
		sleep 0.2
	done
done