mkdir -p crest
cd crest
crest ../traj.xyz --gfnff --alpb chcl3 --quick -nmr --rthr 0.175 --bthr 0.03 --ethr 0.15 | tee  ../crest/traj.out
cp anmr_nucinfo anmr_rotamer coord crest_conformers.xyz ../
cd ..
