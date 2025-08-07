#!/bin/bash
# create a xyz including start and end indices  
# use geodesic_interpolate to generate the interpolated.xyz
# use xyz2NEBallxyz.sh to generate the interpolated.allxyz
# double check the traj.inp is correct or not
# use orca program to run NEB  
#
# reverse operation -  Under vim use this ":g/>/d" to delete any line containing character (>) 

declare -i Rows Nums

Rows=$(sed -n '1p' interpolated.xyz | awk '{print $1+2}')
echo ${Rows}

split -l ${Rows} -d -a 4 interpolated.xyz .part_
echo $">" > .newline

Nums=$(ls .part_* | sort -V | sed  -n '$p' | sed 's/.part_//g'|awk '{print$1+1}')
echo ${Nums}
j=0
Zero_first=$(printf %04d $j)
echo ${Zero_first}
cp .part_${Zero_first} start.xyz
j=$((${Nums}-1))
Zero_end=$(printf %04d $j)
echo ${Zero_end}
cp .part_${Zero_end} final.xyz

rm -rf interpolated.allxyz

for (( i=0; i<${Nums}; i=i+1 ))
	do
	Zero_i=$(printf %04d $i)
    cat .part_${Zero_i} >> interpolated.allxyz
    cat .newline >> interpolated.allxyz
	done

R_Nums=$(awk 'END{print NR}' interpolated.allxyz | awk '{print $1-1}')
echo ${R_Nums}
head -n ${R_Nums} interpolated.allxyz > .tmp && mv .tmp interpolated.allxyz

rm .newline .part_*
echo " Completed"