#!/bin/bash

## Set new directory/simulation type
sim=npt
type=NPT

if [[ ! -e "$sim".gro ]]; then

## Set main directory to current path
gmx_run=$(pwd)
echo "home directory of the gromacs simulations set to $gmx_run"

################################################
########### VOLUME EQUILIBRATION  ##############
################################################

echo "Starting "$type" simulation"

cd $gmx_run

latest_gro=$(ls -t1 *.gro | head -n1)

echo $latest_gro

mkdir "$type"

cp $latest_gro ./"$type"/pre"$sim".gro

cd "$type"

if [[ $1 -eq 0 ]]
	then
    	echo "No arguments supplied"
		nsteps=20000
	else
		nsteps=$1
fi

if [[ $1 -eq 0 ]]
	then
    	echo "No arguments supplied"
		framerate=1000
	else
		framerate=$2
fi


### Create em.mdp file ########################################
cat << EOF >| "$sim".mdp
title		    = $type equilibration 
define      	= -DGMINFF_k500 -DOPC3_IOD_LM
; Run parameters
integrator		= md		; leap-frog integrator
nsteps			= $nsteps	; run these many steps
nstcomm     	= 1000      ; Restrain COM every ... step
dt	        	= 0.001	; timestep, 1 fs
comm-mode       = Linear	; Linear COM mass control every 10 steps
comm-grps       = $commgrps; Groups to apply COM control to 
; Output control
nstxout     	= 0         ; save coordinates (high precision) every 0 timestep
nstvout			= 0		    ; save velocities every 0 timestep
nstfout    		= 0		    ; save forces every 0 timestep
nstxout-compressed		= $framerate		; xtc compressed trajectory output every x steps
nstenergy		= 100		; save energies every 0.1 ps
nstlog			= 100		; update log file every 0.1 ps

; Bond parameters
continuation	= yes		; second dynamics run
;constraint_algorithm = lincs	; holonomic constraints 
constraints		= none  	; none
;lincs_iter		= 1			; accuracy of LINCS
;lincs_order	= 4			; also related to accuracy

; Neighborsearching
ns_type			= grid		; search neighboring grid cells
nstlist			= 20		; 20 fs
rlist			= 1.2		; short-range neighborlist cutoff (in nm)
cutoff-scheme        = Verlet
verlet-buffer-tolerance  = 0.005

; Electrostatics and vdw
coulombtype		= PME		; Particle Mesh Ewald for long-range electrostatics
pme_order		= 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
rcoulomb		= 1.2		; short-range electrostatic cutoff (in nm)
rvdw			= 1.2		; short-range van der Waals cutoff (in nm)

; Temperature coupling is on
tcoupl			= V-rescale	        ; modified Berendsen thermostat
tc-grps         = $tcgrps      	    ; two coupling groups - more accurate
tau_t           = 1 ; 1		        ; time constant, in ps
ref_t           = 298 ;298        	; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl			= Berendsen             ; C-rescale ; Parrinello-Rahman		; Pressure coupling on in NPT
pcoupltype		= anisotropic	        ; Independent pressures in x/y/z directions
nstpcouple      = 20					; default -1 == 100
tau-p			= 2                     ; time constant, in ps
ref-p			= 1.0 1.0 1.0 0.0 0.0 0.0	        ; reference pressure, in bar
compressibility = 4.5E-5 4.5E-5 4.5E-5 0 0 0   ; isothermal compressibility of water, bar^-1
refcoord-scaling = com

; Periodic boundary conditions
pbc				= xyz		        ; 3-D PBC
periodic-molecules = yes            ; Periodic molecules in da house!

; Dispersion correction
DispCorr		= EnerPres	        ; account for cut-off vdW scheme

; energygrps		= $energygrps
	 
EOF

################################################
########### VOLUME EQUILIBRATION  ##############
################################################

$gmx_grompp -f ./"$sim".mdp -c pre"$sim".gro -r pre"$sim".gro -p $gmx_run/$moldir/topol.top -n $gmx_run/$moldir/index.ndx -o "$sim".tpr -pp -po -maxwarn 1

$gmx_mdrun -v -deffnm "$sim"

sleep "$Sleep"s

echo "VOLUME EQUILIBRATION COMPLETE."

if [ -e "$sim".gro ]
		then
  			echo "Volume equilibration may have finished ok!"
  			cp "$sim".gro ../
		fi
 
rm -f \#*

else
		echo "NPT have already simulated"
fi

cd $gmx_run
