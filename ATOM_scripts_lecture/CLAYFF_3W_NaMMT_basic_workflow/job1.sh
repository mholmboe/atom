#!/bin/bash

# Set the cluster specific gromacs command here:
gmx_mdrun="gmx mdrun"
gmx_grompp="gmx grompp"


### Run energy minimization
mkdir EM

cd EM

$gmx_grompp -f ../em.mdp -c ../MMT_clayff/preem.gro -r ../MMT_clayff/preem.gro -p ../MMT_clayff/topol.top -n ../MMT_clayff/index.ndx -o em.tpr -maxwarn 1

$gmx_mdrun -v -deffnm em

cd ..

### Run NVT with position restraints
mkdir NVT

cd NVT

$gmx_grompp -f ../nvt.mdp -c ../EM/em.gro -r ../EM/em.gro -p ../MMT_clayff/topol.top -n ../MMT_clayff/index.ndx -o nvt.tpr -maxwarn 1

$gmx_mdrun -v -deffnm nvt

cd ..

### Run NPT to optimize box size
mkdir NPT

cd NPT

$gmx_grompp -f ../npt.mdp -c ../NVT/nvt.gro -p ../MMT_clayff/topol.top -n ../MMT_clayff/index.ndx -o npt.tpr -maxwarn 1

$gmx_mdrun -v -deffnm npt

cd ..

### Run production run in NVT
mkdir MD

cd MD

$gmx_grompp -f ../md.mdp -c ../NPT/npt.gro -p ../MMT_clayff/topol.top -n ../MMT_clayff/index.ndx -o md.tpr -maxwarn 1

$gmx_mdrun -v -deffnm md

cd ..