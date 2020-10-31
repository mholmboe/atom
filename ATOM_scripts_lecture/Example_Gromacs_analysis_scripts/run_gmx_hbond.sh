#!/bin/bash

# bash script to analyze the hydrogen bonds in the system

type=MD
sim=md

cd $type

echo 0 0 | gmx hbond -f $sim.xtc -s $sim.tpr -n ../SOL_gmx/index.ndx -temp 298 -num hbnum.xvg -hbn hbond.ndx -nhbdist nhbdist.xvg -g hbond.log -dist hbdist.xvg -ang hbang.xvg 

gmx analyze -f hbnum.xvg

xmgrace -nxy hbnum.xvg &

rm -f \#*

