mkdir -p crest
cd crest
crest ../traj.xyz --len 2 --quick -nmr | tee  ../crest/traj.out
cp anmr_nucinfo anmr_rotamer coord crest_conformers.xyz ../
cd ..
