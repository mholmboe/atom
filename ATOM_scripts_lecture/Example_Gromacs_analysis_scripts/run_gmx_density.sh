#!/bin/bash

# bash script to analyze the density of the system

type=MD
sim=md

cd  $type
    
echo 0 | gmx density -f $sim.xtc -s $sim.tpr -o density.xvg \
-d z -sl 200

gmx analyze -f density.xvg  -errbar stddev -aver_start 100

rm -f \#*

xmgrace -nxy density.xvg &

cd ..
