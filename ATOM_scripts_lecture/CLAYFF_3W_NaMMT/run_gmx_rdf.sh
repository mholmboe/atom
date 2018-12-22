#!/bin/bash

# bash script to analyze the radial distribution function between OW – HW
# in blocks using a for loop

type=MD
sim=md

cd  $type

block=200
for N in 1 2 3 4 5;
	do (start=$(($N*$block-$block)) 
		stop=$(($N*$block))
    	echo $start
    	echo $stop       
        echo 2 2 | gmx rdf -f $sim.xtc -s $sim.tpr -n ../SOL_gmx/index.ndx \
        -ref OW -sel HW -bin 0.002 -rmax 1.2 -o rdf_OwHw_$N.xvg \
        -cn rdf_OwHw_cn_$N.xvg -b $start -e $stop
        );
done

rm -f \#*

xmgrace -nxy rdf_*.xvg &

cd ..
