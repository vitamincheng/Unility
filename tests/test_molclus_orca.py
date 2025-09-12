#!/usr/bin/env python
import os
import pytest
import argparse
import censo_ext.molclus_orca as orca
import filecmp
import platform
from pathlib import Path
inFile: Path = Path("tests/data/06.EthylAcetate/02.ORCA_r2SCAN_3C/traj.xyz")
inTemplate_sp_File: Path = Path(
    "tests/data/06.EthylAcetate/02.ORCA_r2SCAN_3C/template_sp.inp")
inTemplate_opt_File: Path = Path(
    "tests/data/06.EthylAcetate/02.ORCA_r2SCAN_3C/template_opt.inp")
outFile: Path = Path("isomers.xyz")
_system = platform.system()


def test_orca_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        orca.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_orca_sp():
    x: dict = {"file": inFile, "template": inTemplate_sp_File,
               "reserve": False, "out": outFile}

    args = argparse.Namespace(**x)
    orca.main(args)

    if _system == "Linux":  # Need 2 min
        compare: Path = Path("tests/compare/orca_isomers_sp.xyz")

    elif _system == "Darwin":  # Need 5 min
        compare: Path = Path("tests/compare/orca_isomers_sp_Darwin.xyz")

    assert filecmp.cmp(args.out, compare)  # type: ignore
    os.remove(args.out)
    import subprocess
    subprocess.call("rm -f 000*.xyz 000*.out 000*.gbw", shell=True)


@pytest.mark.slow
def test_orca_opt():
    x: dict = {"file": inFile, "template": inTemplate_opt_File,
               "reserve": False, "out": outFile}

    args = argparse.Namespace(**x)
    orca.main(args)

    if _system == "Linux":  # Need 2 min
        compare: Path = Path("tests/compare/orca_isomers_opt.xyz")

    elif _system == "Darwin":  # Need 5 min
        compare: Path = Path("tests/compare/orca_isomers_opt_Darwin.xyz")

    assert filecmp.cmp(args.out, compare)  # type: ignore
    os.remove(args.out)
    import subprocess
    subprocess.call("rm -f 000*.xyz 000*.out 000*.gbw", shell=True)


@pytest.mark.slow
def test_orca_opt_default():
    x: dict = {"file": inFile, "template": "template.inp",
               "reserve": False, "out": outFile}

    args = argparse.Namespace(**x)
    orca.main(args)

    if _system == "Linux":  # Need 2 min
        compare: Path = Path("tests/compare/orca_isomers_opt.xyz")

    elif _system == "Darwin":  # Need 5 min
        compare: Path = Path("tests/compare/orca_isomers_opt_Darwin.xyz")

    assert filecmp.cmp(args.out, compare)  # type: ignore
    os.remove(args.out)
    import subprocess
    subprocess.call("rm -f 000*.xyz 000*.out 000*.gbw", shell=True)

# For all marker_fuction or not marker_function
# -v verbose;  -s --capture=no; -m MARKEXPR
#   pytest -v -s tests -m "not slow"
#   pytest -v -s tests -m slow
#
# For specific test function
# --pdb Start the interactive Python debugger on errors
#   pytest -v -s --pdb tests/test_molclus_orca.py::test_orca_opt
#   pytest -v -s --pdb tests/test_molclus_orca.py::test_orca_opt_default
