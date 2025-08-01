#!/usr/bin/env python4
import os
import pytest
import argparse
import censo_ext.molclus_gxtb as gxtb
import filecmp
import platform

inFile = "tests/data/06.EthylAcetate/01.Crest/crest_conformers.xyz"
outFile = "tests/compare/isomers.xyz"
_system = platform.system()


def test_gxtb_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        gxtb.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin", reason="g-xtb only work under Linux")
def test_gxtb_alpb_opt():

    x: dict = {"file": inFile, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": True,
               "method": "gxtb", "out": outFile}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    compare = ""
    if _system == "Darwin":
        compare = "tests/compare/molclus_gxtb_1_Darwin.xyz"
    elif _system == "Linux":
        compare = "tests/compare/molclus_gxtb_1.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


@pytest.mark.skipif(_system == "Darwin", reason="g-xtb only work under Linux")
def test_gxtb_gbsa_opt():
    x: dict = {"file": inFile, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": True,
               "method": "gxtb", "out": outFile}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    compare = ""
    if _system == "Darwin":
        # gxtb is not work in Darwin system
        compare = "tests/compare/molclus_gxtb_2_Darwin.xyz"
    elif _system == "Linux":
        compare = "tests/compare/molclus_gxtb_2.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


@pytest.mark.skipif(_system == "Darwin", reason="g-xtb only work under Linux")
def test_gxtb_alpb():
    x: dict = {"file": inFile, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": False,
               "method": "gxtb", "out": outFile}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    compare = ""
    if _system == "Darwin":
        compare = "tests/compare/molclus_gxtb_3_Darwin.xyz"
    elif _system == "Linux":
        compare = "tests/compare/molclus_gxtb_3.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")

    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


@pytest.mark.skipif(_system == "Darwin", reason="g-xtb only work under Linux")
def test_gxtb_gbsa():
    x: dict = {"file": inFile, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": False,
               "method": "gxtb", "out": outFile}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    compare = ""
    if _system == "Darwin":
        # gxtb is not work in Darwin system
        compare = "tests/compare/molclus_gxtb_4_Darwin.xyz"
    elif _system == "Linux":
        compare = "tests/compare/molclus_gxtb_4.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")

    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)
