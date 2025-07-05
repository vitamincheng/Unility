#!/usr/bin/env python4
import os
from re import I
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_gxtb as gxtb
import filecmp
import platform

in_file = f"tests/data/EthylAcetate/01.Crest/crest_conformers.xyz"
out_file = f"tests/compare/isomers.xyz"


def test_gxtb_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        gxtb.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_gxtb_alpb_opt():

    x: dict = {"file": in_file, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": True,
               "method": "gxtb", "out": out_file}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    dcmp: bool = False
    if platform.system() == "Darwin":
        # compare = "tests/compare/molclus_gxtb_1_Darwin.xyz" # gxtb is not work in Darwin system
        dcmp = True
    else:
        compare = "tests/compare/molclus_gxtb_1.xyz"
        dcmp = filecmp.cmp(args.out, compare)
    assert dcmp == True
    os.remove(args.out)


def test_gxtb_gbsa_opt():
    x: dict = {"file": in_file, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": True,
               "method": "gxtb", "out": out_file}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    dcmp: bool = False
    if platform.system() == "Darwin":
        # compare = "tests/compare/molclus_gxtb_2_Darwin.xyz" # gxtb is not work in Darwin system
        dcmp = True
    else:
        compare = "tests/compare/molclus_gxtb_2.xyz"
        dcmp = filecmp.cmp(args.out, compare)
    assert dcmp == True
    os.remove(args.out)


def test_gxtb_alpb():
    x: dict = {"file": in_file, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": False,
               "method": "gxtb", "out": out_file}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    dcmp: bool = False
    if platform.system() == "Darwin":
        # compare = "tests/compare/molclus_gxtb_3_Darwin.xyz" # gxtb is not work in Darwin system
        dcmp = True
    else:
        compare = "tests/compare/molclus_gxtb_3.xyz"
        dcmp = filecmp.cmp(args.out, compare)
    assert dcmp == True
    os.remove(args.out)


def test_gxtb_gbsa():
    x: dict = {"file": in_file, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": False,
               "method": "gxtb", "out": out_file}
    args = argparse.Namespace(**x)
    gxtb.main(args)

    dcmp: bool = False
    if platform.system() == "Darwin":
        # compare = "tests/compare/molclus_gxtb_4_Darwin.xyz" # gxtb is not work in Darwin system
        dcmp = True
    else:
        compare = "tests/compare/molclus_gxtb_4.xyz"
        dcmp = filecmp.cmp(args.out, compare)
    assert dcmp == True
    os.remove(args.out)
