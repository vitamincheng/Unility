#!/bin/bash
xtb traj.xyz --namespace thermo --ohess
xtb thermo --sthr 100 --temp 298.15 thermo.xtbopt.xyz thermo.hessian > out
rm thermo*