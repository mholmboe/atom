#!/bin/bash

## Set new directory/simulation type
sim=npt
type=NPT

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
define      	= -D$watermodel ;-DPOSRES_noH
; Run parameters
integrator		= md		; leap-frog integrator
nsteps			= $nsteps	; run these many steps
nstcomm     	= 1000      ; Restrain COM every ... step
dt	        	= 0.001	    ; timestep, 1 fs
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
constraint_algorithm = lincs	; holonomic constraints 
constraints		= h-bonds  	; none
lincs_iter		= 1			; accuracy of LINCS
lincs_order	= 4			; also related to accuracy

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
rvdw-switch     = 1.0
vdwtype         = cutoff
vdw-modifier    = force-switch

; Temperature coupling is on
tcoupl			= V-rescale	        ; modified Berendsen thermostat
tc-grps         = $tcgrps      	    ; two coupling groups - more accurate
tau_t           = 2   ;2		    ; time constant, in ps
ref_t           = 298 ;298        	; reference temperature, one for each group, in K

; Pressure coupling is on
pcoupl			= Berendsen             ; Parrinello-Rahman		; Pressure coupling on in NPT
pcoupltype		= anisotropic	        ; keeps same pressure in x/y directions, separate in z
tau-p			= 2            			; time constant, in ps
ref-p			= 1.0 1.0 1.0 0.0 0.0 0.0	        ; reference pressure, in bar
compressibility = 0 4.5E-5 0 0 0 0        ; isothermal compressibility of water, bar^-1
refcoord-scaling = com

; Periodic boundary conditions
pbc				= xyz		; 3-D PBC
periodic-molecules = yes

; Dispersion correction
DispCorr		= no ; EnerPres	; account for cut-off vdW scheme

energygrps		= $energygrps
	 
EOF

################################################
########### VOLUME EQUILIBRATION  ##############
################################################

$run_grompp -f ./"$sim".mdp -c pre"$sim".gro -p $gmx_run/$moldir/topol.top -n $gmx_run/$moldir/index.ndx -o "$sim".tpr -pp -po -maxwarn 1

$run_cmd -v -deffnm "$sim"

sleep "$Sleep"s

echo "VOLUME EQUILIBRATION COMPLETE."

if [ -e "$sim".gro ]
		then
  			echo "Volume equilibration may have finished ok!"
  			cp "$sim".gro ../
		fi
 
rm -f \#*

cd $gmx_run
