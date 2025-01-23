import string
import pytest
from pathlib import Path
import argparse
import xyzTranslate


def test_xyzTranslate_move():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
               "move": [5, 0, 0], "out": "tests/compare/output.xyz", "cut": None}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    import filecmp
    compare = "tests/compare/xyzTranslate-move.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    import os
    os.remove(args.out)


def test_xyzTranslate_cut():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
               "move": [5, 0, 0], "out": "tests/compare/output.xyz", "cut": 10}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    import filecmp
    compare = "tests/compare/xyzTranslate-cut.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    import os
    os.remove(args.out)


def test_xyzTranslate_cut_move():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
               "move": [5, 0, 0], "out": "tests/compare/output.xyz", "cut": 3}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)

    x: dict = {"file": "tests/compare/output.xyz",
               "move": [0, 0, 5], "out": "tests/compare/output2.xyz", "cut": 3}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    import os

    import filecmp
    compare = "tests/compare/xyzTranslate-cut-move.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    import os
    os.remove(args.file)
    os.remove(args.out)
