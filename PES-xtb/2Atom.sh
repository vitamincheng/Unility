#!/bin/bash
mkdir -p xtb_2Atom
cd xtb_2Atom
xtb ../xtbscan_single.xyz --gfn 2 --opt --input ../scan.inp
cd ..
cat xtb_2Atom/xtbscan.log >> xtbscan2.xyz