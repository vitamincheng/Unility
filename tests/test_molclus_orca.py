#!/usr/bin/env python3
import os
import pytest
import argparse
import censo_ext.molclus_orca as orca
import filecmp
import platform

inFile = "tests/data/06.EthylAcetate/02.ORCA_r2SCAN_3C/traj.xyz"
outFile = "isomers.xyz"
_system = platform.system()


def test_orca_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        orca.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin", reason="")
def test_orca_opt():
    x: dict = {"file": inFile, "template": "template.inp",
               "remove": True, "out": outFile}

    args = argparse.Namespace(**x)
    orca.main(args)

    compare = ""
    if _system == "Linux":  # Need 2 min
        compare = "tests/compare/orca_isomers.xyz"

    elif _system == "Darwin":
        compare = "tests/compare/orca_isomers_Darwin.xyz"

    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)
    import subprocess
    subprocess.call("rm -f 000*.xyz 000*.out 000*.gbw", shell=True)
