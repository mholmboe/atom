#!/bin/bash

# Make sure you understand what all the below gromacs commands and -flags do, by for instance running 
# gmx grompp -h or gmx mdrun -h in a separate terminal window

# Set the cluster specific gromacs command here, in case they differ from the default on your computer or HPC cluster.
gmx_grompp="gmx grompp"
gmx_mdrun="gmx mdrun"

### Run ENERGY MINIMIZATION - to relax the system and make sure we do not have overlapping atoms etc.
mkdir EM # First create a folder for energy minimization
cd EM # Change directory to the newly created EM folder

# Run Gromacs pre-processing tool gmx grompp to generate a .tpr run input file.
$gmx_grompp -f ../em.mdp -c ../gmx_minff/preem.gro -r ../gmx_minff/preem.gro -p ../gmx_minff/topol.top -n ../gmx_minff/index.ndx -o em.tpr -pp -po -maxwarn 1  
			
# Run the simulation with Gromacs gmx mdrun, with the general name definition em (as in em.tpr for instance)
$gmx_mdrun -v -deffnm em 

# Reflect over the generated output and the different files created. Did the simulation finish ok?
 
cd .. # Change directory back to the parent folder
###


### Run simulation in NVT with position restraints - mainly to equilibrate the solvent and solute ions
mkdir NVT
cd NVT
$gmx_grompp -f ../nvt.mdp -c ../EM/em.gro -r ../EM/em.gro -p ../gmx_minff/topol.top -n ../gmx_minff/index.ndx -o nvt.tpr -pp -po -maxwarn 1
$gmx_mdrun -v -deffnm nvt
cd ..
###


### Run simulation in NPT - to equilibrate the density of the system by optimizing the box size
mkdir NPT
cd NPT
$gmx_grompp -f ../npt.mdp -c ../NVT/nvt.gro -p ../gmx_minff/topol.top -n ../gmx_minff/index.ndx -o npt.tpr -pp -po -maxwarn 2
$gmx_mdrun -v -deffnm npt
cd ..
###


### Run the PRODUCTION run (here in NVT)
mkdir MD
cd MD
$gmx_grompp -f ../md.mdp -c ../NPT/npt.gro -p ../gmx_minff/topol.top -n ../gmx_minff/index.ndx -o md.tpr -pp -po -maxwarn 1
$gmx_mdrun -v -deffnm md
cd ..
###