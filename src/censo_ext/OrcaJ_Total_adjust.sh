#!/bin/bash
declare -i CONF_Number answer

CONF_Number=$(ls -d CONF* | sed 's/CONF/ /' | sort -g | tail -n1 | awk '{print $1}')

echo ${CONF_Number}
	
for (( i=1; i<=${CONF_Number}; i=i+1 ))
	do
	if [ -d CONF${i}/NMR ]; then
	    sed -i 's/(Total)/ Total /g' CONF${i}/NMR/orcaJ.out
    else echo "CONF${i}/NMR. Directory does not exist." 
	fi
	done
echo " Completed"
