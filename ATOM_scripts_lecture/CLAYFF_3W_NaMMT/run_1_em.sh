#!/bin/bash

## Set new directory/simulation type
sim=em
type=EM

## Set main directory to current path
gmx_run=$(pwd)
echo "home directory of the gromacs simulations set to $gmx_run"

################################################
####### ENERGY MINIMIZATION WITH STEEP #########
################################################

echo "Starting energy minimization steep"

cd $gmx_run

mkdir "$type"

cd "$type"

cp $gmx_run/$moldir/preem.gro ./pre"$sim".gro

if [[ $1 -eq 0 ]]
	then
    	echo "No arguments supplied"
		nsteps=1000
	else
		nsteps=$1
fi

if [[ $1 -eq 0 ]]
	then
    	echo "No arguments supplied"
		framerate=10
	else
		framerate=$2
fi


### Create em.mdp file ########################################
cat << EOF >| "$sim".mdp
; Parameters describing what to do, when to stop and what to save
define		= -DFLEXIBLE -D$watermodel -DPOSRES_noH ; Use felxible SPC/E during miniposition restrain the protein
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 500.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nstxout		= $framerate	; save coordinates every 0 ps
nsteps		= $nsteps 	; Maximum number of (minimization) steps to perform

;freezegrps  = MMT
;freezedim   = Y Y Y

; Bond parameters
constraint_algorithm = lincs	; holonomic constraints 
constraints		= none  	; none
lincs_iter		= 1			; accuracy of LINCS
lincs_order	    = 4			; also related to accuracy

; Neighborsearching
cutoff-scheme        = Verlet

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		= 1 		    ; Frequency to update the neighbor list and long range forces
ns_type		= grid		    ; Method to determine neighbor list (simple, grid)
rlist		= 1.0		    ; Cut-off for making neighbor list (short range forces)
coulombtype	= PME		    ; Treatment of long range electrostatic interactions
rcoulomb	= 1		        ; Short-range electrostatic cut-off
rvdw		= 1		        ; Short-range Van der Waals cut-off
pbc		    = xyz 		    ; Periodic Boundary Conditions
periodic-molecules = yes    ; Periodic molecules in da house!
	 
EOF

################################################
####### ENERGY MINIMIZATION WITH STEEP #########
################################################
echo $gmx_run/$moldir

$gmx_grompp -f ./"$sim".mdp -c pre"$sim".gro -r pre"$sim".gro -p $gmx_run/$moldir/topol.top -n $gmx_run/$moldir/index.ndx -o "$sim".tpr -pp -po -maxwarn 2
 
$gmx_mdrun -v -deffnm "$sim"

sleep "$Sleep"s

echo "Minimization steep complete."

if [ -e "$sim".gro ]
		then
  			echo "Minimization steep may have finished ok!"
  			cp "$sim".gro ../
		fi
 
#rm -f \#*

cd $gmx_run