#!/bin/bash

# bash script to analyze the density of the system

type=MD
sim=md

cd  $type

echo 0 | gmx msd -f md.xtc -s md.tpr -o msd.xvg -trestart 250 -beginfit 20 -endfit 250

xmgrace -nxy msd.xvg &


