#!/bin/sh

#PBS -M borisl@technion.ac.il

#PBS -mbea

#PBS -l select=1:ncpus=1 -l place=scatter

PBS_O_WORKDIR=$HOME/Iliya/microgrids
cd $PBS_O_WORKDIR
now=$(date +"../logs/%Y-%m-%d_%H-%M-%S.%3N_output.log")

#MATLAB run command
#-----------------------
#Assuming the script is called using: qsub -v fid=5,iter=404,psf=0.75 run_sim.sh
matlab2020a -batch "runAlgorithms_1p(${fid}, ${iter}, ${psf})" > ${now}


