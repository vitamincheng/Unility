#!/usr/bin/env python3
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzSplit as xyzSplit
import filecmp
import os


def test_xyzSplit():
    x: dict = {"file": "tests/data/crest_conformers1.xyz", "cuts": 12,
               "atoms": [52, 55], "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)
    xyzSplit.main(args)

    compare = "tests/compare/xyzSplit.xyz"
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)


def test_xyzSplit_miss_cuts():
    x: dict = {"file": "tests/data/crest_conformers1.xyz", "cuts": None,
               "atoms": [52, 55], "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(args)
    assert e.type == SystemExit
    assert e.value.code == 0    # for argprarse wrong


def test_xyzSplit_miss_atoms():
    x: dict = {"file": "tests/data/crest_conformers1.xyz", "cuts": 12,
               "atoms": None, "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(args)
    assert e.type == SystemExit
    assert e.value.code == 0    # for argprarse wrong


def test_xyzSplit_miss_cuts_atoms():
    x: dict = {"file": "tests/data/crest_conformers1.xyz", "cuts": None,
               "atoms": None, "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(args)
    assert e.type == SystemExit
    assert e.value.code == 0    # for argprarse wrong


def test_xyzSplit_miss():
    x: dict = {}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        xyzSplit.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2    # for argprarse wrong
