import string
import pytest
from pathlib import Path
import argparse
import xyzSplit


def test_xyzSplit():
    x: dict = {"file": "tests/data/crest_conformers1.xyz", "cut": 12,
               "atom": [52, 55], "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)
    xyzSplit.main(args)

    import filecmp
    compare = "tests/compare/xyzSplit.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    import os
    os.remove(args.out)
