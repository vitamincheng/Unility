#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_tblite_xtb as tblite_xtb
import filecmp
import platform

inFile = f"tests/data/06.EthylAcetate/01.Crest/crest_conformers.xyz"
outFile = f"tests/compare/isomers.xyz"
_system = platform.system()


def test_tblite_xtb_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        tblite_xtb.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin" or _system == "Linux", reason="tblite xtb have some bugs under Darwin and Linux")
def test_tblite_xtb_alpb_opt():

    x: dict = {"file": inFile, "alpb": "chloroform", "gbsa": None, "chrg": 0, "uhf": 1, "opt": True,
               "method": "GFN2-xTB", "out": outFile}
    args = argparse.Namespace(**x)
    tblite_xtb.main(args)

    dcmp: bool = False
    compare = ""
    if _system == "Darwin":
        compare = f"tests/compare/molclus_xtb_1_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_tblite_xtb_1.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    dcmp = filecmp.cmp(args.out, compare)
    # os.remove(args.out)
    assert dcmp == True


@pytest.mark.skipif(_system == "Darwin" or _system == "Linux", reason="tblite xtb have some bugs under Darwin and Linux")
def test_tblite_xtb_gbsa_opt():
    x: dict = {"file": inFile, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2", "out": outFile}
    args = argparse.Namespace(**x)
    tblite_xtb.main(args)

    dcmp: bool = False
    compare = ""
    if _system == "Darwin":
        compare = f"tests/compare/molclus_xtb_2_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_xtb_2.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    dcmp = filecmp.cmp(args.out, compare)
    os.remove(args.out)
    assert dcmp == True


def test_tblite_xtb_alpb():
    x: dict = {"file": inFile, "alpb": "chloroform", "gbsa": None, "chrg": 0, "uhf": 1, "opt": False,
               "method": "GFN2-xTB", "out": outFile}
    args = argparse.Namespace(**x)
    tblite_xtb.main(args)

    dcmp = False
    compare = ""
    if _system == "Darwin":
        compare = f"tests/compare/molclus_tblite_xtb_3_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_tblite_xtb_3.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    dcmp = filecmp.cmp(args.out, compare)
    os.remove(args.out)
    assert dcmp == True


@pytest.mark.skipif(_system == "Darwin" or _system == "Linux", reason="tblite xtb have some bugs under Darwin and Linux")
def test_tblite_xtb_gbsa():
    x: dict = {"file": inFile, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": False,
               "method": "gfn2", "out": outFile}
    args = argparse.Namespace(**x)
    tblite_xtb.main(args)

    dcmp = False
    compare = ""
    if _system == "Darwin":
        compare = f"tests/compare/molclus_xtb_4_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_xtb_4.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    dcmp = filecmp.cmp(args.out, compare)
    os.remove(args.out)
    assert dcmp == True
