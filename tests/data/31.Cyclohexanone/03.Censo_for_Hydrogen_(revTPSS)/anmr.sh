#!/bin/bash 			
#		                 	[06.24.2023] vitamin.cheng@gmail.com
function anmrh () {
#	cp -rf .anmrrcH .anmrrc
	anmr -mss 10 | tee anmrh.out
	head -1 anmr.dat > ~anmr.dat
	awk '($2 > 0.001){print $0}' anmr.dat >> ~anmr.dat
	tail -1 anmr.dat >> ~anmr.dat
	mv -f ~anmr.dat anmrh.dat
	rm anmr.dat
}

function anmrhcompact () {
	
	awk '{printf "    %7.4f    %12.6e\n",$1,$2}' anmrh.dat |  awk '!visited[$1]++ { print $0 }' > ~anmrh.dat
	mv -f ~anmrh.dat anmrh.dat
}

function anmrhmagnitude () {
	
	magnitude=100000000

	if [ ! -z $1 ] && [ -f $1 ] ; then
		magnitude=$(sort -g -rk2,2 $1 | head -1 | awk '{print $2}' )
	else echo "The default of magnitude is 1.00+e8"
	fi	
	echo $magnitude
	anmrh_magnitude=$(sort -g -rk2,2 anmrh.dat | head -1 | awk '{print $2}' )
	echo $anmrh_magnitude
	rescale=$(echo "${anmrh_magnitude} ${magnitude}" | awk '{print $2/$1}')
	echo $rescale
	
	awk '{printf "    %7.4f    %12.6e\n",$1,$2*v1}' v1=$rescale  anmrh.dat > ~anmrh.dat
	
	mv -f ~anmrh.dat anmrh.dat
}



function Help () {
	echo "================================================================================" 
	echo "* 					   [01.06.2023] vitamin.cheng@gmail.com"
	echo "* For anmr from censo program (Modification)(only for revTPSS)                  "
	echo "* Input  :                                                                     "
	echo "* Needed : .anmrrc anmr_enso anmr_nucinfo (from crest) 		               "
	echo "* Output : anmrh.dat andd anmrh.out    "
	echo "* Option :    "
	echo "*    -h          - only for proton spectra "
	echo "*    -hcompact   - reduce the size of anmrh.dat"
	echo "*    -hmagnitude - rescale magnitude to maxium height of reference" 
	echo "*                  (default reference is 100000000)"
	echo "================================================================================"
}



if [ $# == 0 ]; then
	anmrh
else
	case "$1" in
		-h) echo "Outfile is anmrh.dat and anmrh.out"
			anmrh
			;;
		-hcompact) echo "Outfile is anmrh.dat(compact)"
			anmrhcompact
			;;
		-hmagnitude) echo "Outfile is anmrh.dat"
			anmrhmagnitude "$2"
			;;
                -hsubstraction) echo "Outfile is diff.dat"
                        anmrhsubstraction "$2"
                        ;;
		--help) Help
			;;
		*) echo "only for h and c spectrum -h, -hcompact, -c, -hc and -x"
			"More information : --help"
	esac
fi

#p -rf .anmrrcC .anmrrc
#./anmr -lw 20 | tee anmrc.out
#mv anmr.dat anmrc.dat


#./anmr | tee anmr.out
