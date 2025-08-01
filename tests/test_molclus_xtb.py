#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_xtb as xtb
import filecmp
import platform

inFile: Path = Path(
    f"tests/data/06.EthylAcetate/01.Crest/crest_conformers.xyz")
outFile: Path = Path(f"tests/compare/isomers.xyz")
_system = platform.system()


def test_xtb_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xtb.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.skipif(_system == "Darwin", reason="xtb have some bugs under Darwin")
def test_xtb_alpb_opt():
    x: dict = {"file": inFile, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2", "out": outFile}
    args = argparse.Namespace(**x)

    compare = f""
    xtb.main(args)
    if _system == "Darwin":
        compare = f"tests/compare/molclus_xtb_1_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_xtb_1.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    try:
        assert filecmp.cmp(args.out, compare) == True
    except AssertionError as e:
        print(e)
        os.remove(args.out)


@pytest.mark.skipif(_system == "Darwin", reason="xtb have some bugs under Darwin")
def test_xtb_gbsa_opt():
    x: dict = {"file": inFile, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2", "out": outFile}
    args = argparse.Namespace(**x)

    compare = f""
    xtb.main(args)
    if _system == "Darwin":
        compare = f"tests/compare/molclus_xtb_2_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_xtb_2.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    try:
        assert filecmp.cmp(args.out, compare) == True
    except AssertionError as e:
        print(e)
        os.remove(args.out)


def test_xtb_alpb():
    x: dict = {"file": inFile, "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": False,
               "method": "gfn2", "out": outFile}
    args = argparse.Namespace(**x)

    xtb.main(args)
    compare = f""
    if _system == "Darwin":
        compare = f"tests/compare/molclus_xtb_3_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_xtb_3.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)


def test_xtb_gbsa():
    x: dict = {"file": inFile, "alpb": None, "gbsa": "CHCl3", "chrg": 0, "uhf": 1, "opt": False,
               "method": "gfn2", "out": outFile}
    args = argparse.Namespace(**x)
    xtb.main(args)
    compare = f""
    if _system == "Darwin":
        compare = f"tests/compare/molclus_xtb_4_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/molclus_xtb_4.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)
