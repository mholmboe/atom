#!/bin/bash

cluster=$1		# Passed variable 1 (machine) from caller, leave empty if local machine
ncores=$2		# Passed variable 2 (number of proc.) from caller, leave empty if local machine

## RUN THE SIMULATION SCRIPTS ####################################
source ./run_0_settings.sh $cluster $ncores
source ./run_1_em.sh 10000 100            # Run energy minimization
source ./run_2_nvt.sh 5000 1000           # Run position restraints in an NVT simulation, i.e. we run at constant n_particles, volume and temperature
source ./run_3_npt.sh 50000 10000         # Run NPT simulation, i.e. we optimize volume by running at constant pressure 
source ./run_4_imd.sh 100000 100           # Run the MD equilibration simulation, i.e. we run NVT at constant n_particles, volume and temperature

