#!/usr/bin/env python
import pytest
import argparse
import censo_ext.xyzSplit as xyzSplit
import filecmp
import os
from pathlib import Path
inFile: Path = Path("tests/data/crest_conformers1.xyz")
outFile: Path = Path("tests/compare/output.xyz")


def test_xyzSplit() -> None:
    x: dict = {"file": inFile, "cuts": 12, "check": False,
               "atoms": [52, 55], "out": outFile, "print": False}
    args = argparse.Namespace(**x)
    xyzSplit.main(args)

    compare: Path = Path("tests/compare/xyzSplit.xyz")
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


def test_xyzSplit_miss_cuts() -> None:
    x: dict = {"file": inFile, "cuts": None, "check": True,
               "atoms": [52, 55], "out": outFile, "print": False}

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 0    # for argprarse wrong


def test_xyzSplit_miss_atoms() -> None:
    x: dict = {"file": inFile, "cuts": 12, "check": True,
               "atoms": None, "out": outFile, "print": False}

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 0    # for argprarse wrong


def test_xyzSplit_miss_cuts_atoms() -> None:
    x: dict = {"file": inFile, "cuts": None, "check": True,
               "atoms": None, "out": outFile, "print": False}

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 0    # for argprarse wrong


def test_xyzSplit_miss() -> None:
    x: dict = {}

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2    # for argprarse wrong
