#!/bin/bash 			
#		                 	[07.26.2025] vitamin.cheng@gmail.com
set -e
OrcaJ_Total_adjust.sh
anmr -mss 10 | tee anmr.out
head -1 anmr.dat > ~anmr.dat
awk '($2 > 0.001){print $0}' anmr.dat >> ~anmr.dat
tail -1 anmr.dat >> ~anmr.dat
rm anmr.dat
mv -f ~anmr.dat anmr.dat
rm tmpanmr*