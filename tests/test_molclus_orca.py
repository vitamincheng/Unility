#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_orca as orca
import filecmp
import platform

in_file = f"tests/data/EthylAcetate/02.ORCA_r2SCAN_3C/traj.xyz"
out_file = f"isomers.xyz"
_system = platform.system()


def test_orca_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        orca.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin", reason="")
def test_orca_opt():
    x: dict = {"file": in_file, "template": "template.inp",
               "remove": True, "out": out_file}

    args = argparse.Namespace(**x)
    orca.main(args)

    compare = ""
    if _system == "Linux":  # Need 2 min
        compare = f"tests/compare/orca_isomers.xyz"

    elif _system == "Darwin":
        compare = f"tests/compare/orca_isomers_Darwin.xyz"

    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)
    import subprocess
    subprocess.call("rm -f 000*.xyz 000*.out 000*.gbw", shell=True)
