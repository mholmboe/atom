#!/bin/bash

## Make the index file, in case you forgot #######################
echo q | gmx make_ndx -f preem.gro

# cat << EOF | gmx make_ndx -f preem.gro
# 
# r 1
# 
# name 7 MMT_1
# 
# r 2
# 
# name 8 MMT_2
# 
# 2 & a H
# 
# name 9 Ho
# 
# a Ow* OW*
# 
# name 10 Ow
# 
# a Hw* HW*
# 
# name 11 Hw
# 
# 0 &! 2
# 
# name 12 non-MMT
# 
# q
# 
# EOF

cat << EOF | gmx make_ndx -f preem.gro

r 1

name 9 MMT_1

r 2

name 10 MMT_2

2 & a H

name 11 Ho

a Ow* OW*

name 12 Ow

a Hw* HW*

name 13 Hw

0 &! 2

name 14 non-MMT

q

EOF

rm -f \#*
