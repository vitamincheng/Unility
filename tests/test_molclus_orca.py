#!/usr/bin/env python3
import os
from re import sub
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_orca as orca
import filecmp
import platform
in_file = f"tests/data/EthylAcetate/02.ORCA_r2SCAN_3C/traj.xyz"
out_file = f"isomers.xyz"


def test_orca_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        orca.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_orca_opt():
    x: dict = {"file": in_file, "template": "template.inp",
               "remove": True, "out": out_file}

    args = argparse.Namespace(**x)
    if platform.system() != "Darwin":
        orca.main(args)

    compare = f"tests/compare/orca_isomers.xyz"
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)
    import subprocess
    subprocess.call("rm -f 000*.xyz 000*.out 000*.gbw", shell=True)
