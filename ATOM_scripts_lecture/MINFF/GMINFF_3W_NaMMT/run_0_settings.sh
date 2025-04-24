#!/bin/bash

cluster=$1		# Passed variable 1 (machine) from caller
ncores=$2		# Passed variable 2 (number of proc.) from caller

echo $nPME

if [[ $cluster == "somecluster" ]]
	then
		echo "Cluster is cori"
        echo "will run on $ncores proc."
	
gmx_mdrun="srun -n ${ncores} mdrun_mpi_sp -tunepme" #-npme ${nPME}"

gmx_tune_pme="gmx_sp tune_pme"
export MDRUN='mdrun_mpi_sp'
export MPIRUN='srun'

gmx_convert_tpr="gmx_sp convert-tpr"
gmx_density="gmx_sp density"
gmx_distance="gmx_sp distance"
gmx_distance="gmx_sp distance"
gmx_editconf="gmx_sp editconf"
gmx_grompp="gmx_sp grompp"
gmx_hbond="gmx_sp hbond"
gmx_make_ndx="gmx_sp make_ndx"
gmx_msd="gmx_sp msd"
gmx_rdf="gmx_sp rdf"
gmx_trjconv="gmx_sp trjconv"
gmx_wham="gmx_sp wham"
        
else 

echo "no cluster, running on local machine"
echo "will run on $ncores proc."
	
gmx_mdrun="gmx mdrun"
	
gmx_convert_tpr="gmx convert-tpr"
gmx_density="gmx density"
gmx_distance="gmx distance"
gmx_distance="gmx distance"
gmx_editconf="gmx editconf"
gmx_grompp="gmx grompp"
gmx_hbond="gmx hbond"
gmx_make_ndx="gmx make_ndx"
gmx_msd="gmx msd"
gmx_rdf="gmx rdf"
gmx_trjconv="gmx trjconv"
gmx_wham="gmx wham"
    
fi

echo $gmx_mdrun
echo $gmx_grompp
