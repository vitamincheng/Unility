#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_xtb as xtb
import filecmp
import platform

in_file = f"tests/data/EthylAcetate/01.Crest/crest_conformers.xyz"
out_file = f"tests/compare/isomers.xyz"


def test_xtb_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xtb.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_xtb_alpb_opt():

    x: dict = {"file": in_file, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2", "out": out_file}
    args = argparse.Namespace(**x)

    dcmp: bool = False
    if platform.system() == "Darwin":
        # xtb.main(args) # Under Darwin system fortrans of xtb have some bugs, only for opt
        # compare = "tests/compare/molclus_xtb_1_Darwin.xyz"
        dcmp = True
    elif platform.system() == "Linux":
        xtb.main(args)
        compare = "tests/compare/molclus_xtb_1.xyz"
        dcmp = filecmp.cmp(args.out, compare)
        os.remove(args.out)
    else:
        pass
    assert dcmp == True


def test_xtb_gbsa_opt():
    x: dict = {"file": in_file, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2", "out": out_file}
    args = argparse.Namespace(**x)

    dcmp: bool = False
    if platform.system() == "Darwin":
        # xtb.main(args) # Under Darwin system fortrans of xtb have some bugs, only for opt
        # compare = "tests/compare/molclus_xtb_2_Darwin.xyz"
        dcmp = True
    elif platform.system() == "Linux":
        xtb.main(args)
        compare = "tests/compare/molclus_xtb_2.xyz"
        dcmp = filecmp.cmp(args.out, compare)
        os.remove(args.out)
    else:
        pass
    assert dcmp == True


def test_xtb_alpb():
    x: dict = {"file": in_file, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": False,
               "method": "gfn2", "out": out_file}
    args = argparse.Namespace(**x)
    xtb.main(args)

    if platform.system() == "Darwin":
        compare = "tests/compare/molclus_xtb_3_Darwin.xyz"
    elif platform.system() == "Linux":
        compare = "tests/compare/molclus_xtb_3.xyz"
    else:
        compare = ""
    dcmp = filecmp.cmp(args.out, compare)
    assert dcmp == True
    os.remove(args.out)


def test_xtb_gbsa():
    x: dict = {"file": in_file, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": False,
               "method": "gfn2", "out": out_file}
    args = argparse.Namespace(**x)
    xtb.main(args)

    if platform.system() == "Darwin":
        compare = "tests/compare/molclus_xtb_4_Darwin.xyz"
    elif platform.system() == "Linux":
        compare = "tests/compare/molclus_xtb_4.xyz"
    else:
        compare = ""
    dcmp = filecmp.cmp(args.out, compare)
    assert dcmp == True
    os.remove(args.out)
