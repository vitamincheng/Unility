#!/bin/bash
mkdir -p xtb_1Atom
cd xtb_1Atom
xtb ../traj.xyz --gfn 2 --opt --input ../scan.inp
cd ..
cp xtb_1Atom/xtbscan.log xtbscan.xyz
rm -rf xtbscan2.xyz