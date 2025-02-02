import os
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_xtb as xtb
import filecmp
import platform


def test_xtb_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xtb.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_xtb_alpb_opt():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2", "out": "tests/compare/isomers.xyz"}
    args = argparse.Namespace(**x)
    xtb.main(args)

    compare = "tests/compare/molclus_xtb_1.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    if platform.system() != 'Darwin':
        assert (dcmp == True)
    os.remove(args.out)


def test_xtb_gbsa_opt():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2", "out": "tests/compare/isomers.xyz"}
    args = argparse.Namespace(**x)
    xtb.main(args)

    compare = "tests/compare/molclus_xtb_2.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    if platform.system() != 'Darwin':
        assert (dcmp == True)

    os.remove(args.out)


def test_xtb_alpb():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": False,
               "method": "gfn2", "out": "tests/compare/isomers.xyz"}
    args = argparse.Namespace(**x)
    xtb.main(args)

    compare = "tests/compare/molclus_xtb_3.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    if platform.system() != 'Darwin':
        assert (dcmp == True)
    os.remove(args.out)


def test_xtb_gbsa():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": False,
               "method": "gfn2", "out": "tests/compare/isomers.xyz"}
    args = argparse.Namespace(**x)
    xtb.main(args)

    compare = "tests/compare/molclus_xtb_4.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    if platform.system() != 'Darwin':
        assert (dcmp == True)
    os.remove(args.out)
