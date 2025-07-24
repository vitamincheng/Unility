#!/bin/bash 				
#									[06.24.2023] vitamin.cheng@gmail.com
function anmrh () {
    rm -rf tmp
    anmr -mss 10 -ascal 0.91348 -bscal 28.7564 | tee anmrh.out
	
    head -1 anmr.dat > ~anmr.dat
	awk '($2 > 0.001){print $0}' anmr.dat >> ~anmr.dat
	tail -1 anmr.dat >> ~anmr.dat
	mv -f ~anmr.dat anmrh.dat
	
    rm -rf anmr.dat
}


function anmrc () {

    anmr -mss 10 | tee anmrc.out
    head -1 anmr.dat > ~anmr.dat
    awk '($2 > 0.001){print $0}' anmr.dat >> ~anmr.dat
    tail -1 anmr.dat >> ~anmr.dat
    mv -f ~anmr.dat anmrc.dat
    rm anmr.dat

#    rm -rf tmp
#	anmr -mss 10 -ascal 0.98356 -bscal 195.987 | tee anmrc.out
#	
#	head -1 anmr.dat > ~anmr.dat
#	awk '($2 > 0.000001){print $0}' anmr.dat >> ~anmr.dat
#	tail -1 anmr.dat >> ~anmr.dat
#	mv -f ~anmr.dat anmrc.dat
#    
#    rm -rf anmr.dat
}

function anmrx () {
    rm -rf tmp
    anmr -mss 10 -ascal 0.91348 -bscal 28.7564 | tee anmrx.out
    rm -rf anmr.dat
}


function Help () {
	echo "================================================================================" 
	echo "* 					   [06.21.2023] vitamin.cheng@gmail.com"
	echo "* For anmr from censo program (Modification)(only for PBE0)	            "
	echo "* Input  :                                                                    "
	echo "* Needed : .anmrrc anmr_enso anmr_nucinfo (from crest) 		            "
	echo "* Output : anmrh.dat/anmrc.dat and anmrh.out/anmrc.out/anmrx.out    "
	echo "================================================================================"
}

if [ $# == 0 ]; then
	anmrc
else
	case "$1" in
		-h) echo "Outfile is anmrh.dat and anmrh.out"
			anmrh
			;;
		-c) echo "Outfile is anmrc.dat and anmrc.out"
			anmrc
			;;
		-x) echo "Outfile is anmrx.dat and anmrx.out" 
			anmrh
			anmrc
            anmrx
			;;
		--help) Help
			;;
		*) echo "only for h and c spectrum -h, -c and -x"
			"More information : --help"
	esac
fi

#p -rf .anmrrcC .anmrrc
#./anmr -lw 20 | tee anmrc.out
#mv anmr.dat anmrc.dat


#./anmr | tee anmr.out
