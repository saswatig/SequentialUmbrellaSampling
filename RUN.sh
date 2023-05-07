#!/bin/bash

mkdir RESULTS
rm Codes_UMBRELLA/OUTPUTFILES
	cp Codes_UMBRELLA/INITIAL_PARAMETERS.inc RESULTS/INITIAL_PARAMETERS
	echo ../RESULTS/IConfig>>Codes_UMBRELLA/OUTPUTFILES
	echo ../RESULTS/Config>>Codes_UMBRELLA/OUTPUTFILES
	echo ../RESULTS/Icycle_Energy_Pressure>>Codes_UMBRELLA/OUTPUTFILES
	echo ../RESULTS/FConfig>>Codes_UMBRELLA/OUTPUTFILES
	echo ../RESULTS/Histo_UM>>Codes_UMBRELLA/OUTPUTFILES
cd Codes_UMBRELLA
ifort -O3 *.f -o MC_000.out
#gfortran *.f -o MC.out
./MC_000.out
cd ../
cd POST_PROCESS
bash Script.sh
cd NORMALIZE
bash Script_Norm.sh

