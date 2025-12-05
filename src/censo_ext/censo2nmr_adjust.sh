#!/bin/bash
# censo 2 version
# normal condtion the censo will generate the nmr_s.out and nmr_j.out files
# these files do not read by anmr
# this program will change the filename

declare -i CONF_Number answer

CONF_Number=$(ls -d CONF* | sed 's/CONF/ /' | sort -g | tail -n1 | awk '{print $1}')

echo ${CONF_Number}

for (( i=1; i<=${CONF_Number}; i=i+1 ))
	do
	if [ -f CONF${i}/nmr/nmr_s.out ]; then
	    mv CONF${i}/nmr/nmr_s.out CONF${i}/NMR/orcaS.out
    else echo "CONF${i}/nmr/nmr_s.out The file does not exist." 
	fi
    if [ -f CONF${i}/nmr/nmr_j.out ]; then
	    mv CONF${i}/nmr/nmr_j.out CONF${i}/NMR/orcaJ.out
    else echo "CONF${i}/nmr/nmr_j.out The file does not exist." 
	fi
	done
echo " Completed"
