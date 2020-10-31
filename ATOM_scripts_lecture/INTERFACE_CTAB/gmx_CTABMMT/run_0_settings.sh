#!/bin/bash

## SET TWO MOLECULE/FILENAME KEYWORDS ############################
I=MMT # Keyword describing system
J=interface # Keyword describing forcefield
moldir="$I"_"$J"

Sleep=0		# sleep/take a break in seconds between gromacs runs in case the cluster is slow and stupid.

## SET MAIN DIRECTORY TO CURRENT PATH ############################
gmx_run=$(pwd)
echo "home directory of the gromacs simulations set to $gmx_run"

## MDP SETTINGs ##################################################
energygrps="MMT Water Ion CTA"
tcgrps="System"
commgrps="MMT_1 MMT_2 non-MMT"
watermodel="TIP3P"

cluster=$1		# Passed variable 1 (machine) from caller
ncores=$2		# Passed variable 2 (number of proc.) from caller

if [[ $# -eq 0 ]]; then
   echo "no passed variables"
else
   nPME=$((${ncores}/4)) # Number of dedicated PME cores
fi

echo $nPME

if [[ $cluster == "beskow" ]] # This now seems to be for gromacs 5.1.2. 
    then
        echo "Cluster is beskow"
        echo "will run on $ncores proc."
        
        export OMP_NUM_THREADS=1
        ## APRUN_OPTIONS="-n ${ncores} -d 1 -cc none"
        run_cmd="aprun -n ${ncores} -d 1 -cc none gmx_mpi mdrun -npme ${nPME}"
        run_grompp="gmx_seq grompp"
        run_g_wham="gmx_seq wham"
        run_g_dist="gmx_seq distance"
        run_trjconv="gmx_seq trjconv"
        run_editconf="gmx_seq editconf"
        run_make_ndx="gmx_seq make_ndx"
        run_tpbconv="gmx_seq convert-tpr"

elif [ $cluster = "abisko" ] # This now seems to be for gromacs 5.1.1. 
	then
		echo "Cluster is abisko"
		echo "will run on $ncores proc."
		# For Abisko only! Automatic selection of options to mdrun_mpi depending on parameters given to the SBATCH run command
		if [ -n "$SLURM_CPUS_PER_TASK" ]; then
    		md="mdrun_mpi -ntomp $SLURM_CPUS_PER_TASK -npme ${nPME}"
		else
    		md="mdrun_mpi -ntomp 1"
		fi
        run_cmd="srun $md"
		run_grompp="gmx grompp"
		run_g_wham="gmx wham"
		run_g_dist="gmx distance"	
		run_trjconv="gmx trjconv"
		run_editconf="gmx editconf"
		run_make_ndx="gmx make_ndx"
		run_tpbconv="gmx convert-tpr"
		
elif [[ $cluster = "edison" ]]
	then
		echo "Cluster is edison"
		echo "will run on $ncores proc."
                
        run_cmd="srun -n ${ncores} mdrun_mpi_sp -npme ${nPME}"
        run_grompp="gmx_sp grompp"
        run_g_wham="gmx_sp wham"
        run_g_dist="gmx_sp distance"
        run_trjconv="gmx_sp trjconv"
        run_tpbconv="gmx_sp convert-tpr"

elif [[ $cluster = "cori" ]]
	then
		echo "Cluster is cori"
        echo "will run on $ncores proc."
	
		run_cmd="srun -n ${ncores} mdrun_mpi_sp -npme ${nPME}"
		run_grompp="gmx_sp grompp"
		run_g_wham="gmx_sp wham"
		run_g_dist="gmx_sp distance"	
		run_trjconv="gmx_sp trjconv"
        run_tpbconv="gmx_sp convert-tpr"
    
else 
	echo "no cluster, running on local machine"
	echo "will run on $ncores proc."
	
	# source /usr/local/gromacs-5.1.1/bin/GMXRC
	
	run_cmd="gmx mdrun"
	run_grompp="gmx grompp"
	run_g_wham="gmx wham"
	run_g_dist=" gmx distance"
	run_trjconv="gmx trjconv"
    run_tpbconv="gmx tpbconv"
fi

echo $run_cmd
echo $run_grompp