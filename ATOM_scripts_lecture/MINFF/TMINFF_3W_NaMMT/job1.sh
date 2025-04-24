#!/bin/bash

## SET TWO MOLECULE/FILENAME KEYWORDS ############################
I=gmx # Keyword describing system
J=minff # Keyword describing forcefield
moldir="$I"_"$J"

## SET MAIN DIRECTORY TO CURRENT PATH ############################
gmx_run=$(pwd)
echo "home directory of the gromacs simulations set to $gmx_run"

## MDP SETTINGs ##################################################
energygrps="MMT Ion Water"
tcgrps="System"
commgrps="MMT_1 MMT_2 Water_and_ions"

cluster=$1		# Passed variable 1 (machine) from caller, leave empty if local machine
ncores=$2		# Passed variable 2 (number of proc.) from caller, leave empty if local machine

## RUN THE SIMULATION SCRIPTS ####################################
source ./run_0_settings.sh $cluster $ncores
source ./run_1_em.sh 1000 100              # Run energy minimization
source ./run_2_nvt.sh 50000 1000             # Run position restraints in an NVT simulation, i.e. we run at constant n_particles, volume and temperature
source ./run_3_npt.sh 100000 1000          # Run NPT simulation, i.e. we optimize volume by running at constant pressure 
source ./run_4_md.sh 1000000 1000         # Run the MD equilibration simulation, i.e. we run NVT at constant n_particles, volume and temperature
 
