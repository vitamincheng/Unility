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
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    os.remove(args.out)
