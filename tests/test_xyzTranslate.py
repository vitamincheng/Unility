#!/usr/bin/env python3
import string
from _pytest import monkeypatch
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzTranslate as xyzTranslate
import filecmp
import os


def test_xyzTranslate_move():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
               "move": [5, 0, 0], "out": "tests/compare/output.xyz", "cut": None}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    compare = "tests/compare/xyzTranslate-move.xyz"
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)


def test_xyzTranslate_cut():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
               "move": [5, 0, 0], "out": "tests/compare/output.xyz", "cut": 10}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    compare = "tests/compare/xyzTranslate-cut.xyz"
    assert filecmp.cmp(args.out, compare) == True
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

    compare = "tests/compare/xyzTranslate-cut-move.xyz"
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.file)
    os.remove(args.out)


def test_xyzTranslate_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_xyzTranslate_without_move():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
               "out": "tests/compare/output.xyz", "cut": 10}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(args)
    assert e.type == SystemExit
    assert e.value.code == 1


def test_xyzTranslate_miss_file():
    x: dict = {"file": "tests/data/crest_conformers000.xyz",
               "out": "tests/compare/output.xyz", "cut": 10}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(args)
    assert e.type == SystemExit
    assert e.value.code == 1


# if __name__ == "__main__":
#    import cProfile
#    cProfile.run("test_xyzTranslate_move()", sort="cumtime")
