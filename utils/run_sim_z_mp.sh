#!/bin/sh

#PBS -M borisl@technion.ac.il

#PBS -mbea

PBS_O_WORKDIR=$HOME/Iliya/microgrids
cd $PBS_O_WORKDIR
now=$(date +"../logs/%Y-%m-%d_%H-%M-%S.%3N_output.log")

#MATLAB run command
#-----------------------
#Assuming the script is called using: qsub -v fid=5,minid=500,maxid=510,W=12,psf=0.75 run_sim_z_mp.sh
matlab2020a -batch "warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary'); runAlgorithms(${fid}, ${minid}, ${maxid}, ${W}, ${psf})" > ${now}
