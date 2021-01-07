#!/bin/bash

## Make the index file, in case you forgot #######################
echo q | gmx make_ndx -f preem.gro

cat << EOF | gmx make_ndx -f preem.gro

r 1

name 11 MMT_1

r 2

name 12 MMT_2

q

EOF

# cat << EOF | gmx make_ndx -f preem.gro
# 
# a Na* | a Ca* | a Cl*
# 
# name 3 Ion
# 
# a Ow* OW* Hw* HW* 
# 
# name 4 Water
# 
# 3 | 4
# 
# name 5 Water_and_ions
# 
# a Na
# 
# name 6 Na
# 
# 0 & ! 5
# 
# name 7 MMT
# 
# a H
# 
# name 8 Ho
# 
# a Ow* OW*
# 
# name 9 Ow
# 
# a Hw* HW*
# 
# name 10 Hw
# 
# del 2
# 
# q
# 
# EOF

rm -f \#*
