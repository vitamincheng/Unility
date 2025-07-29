#!/usr/bin/env python3
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzTranslate as xyzTranslate
import filecmp
import os
in_file = f"tests/data/crest_conformers.xyz"
out_file = f"tests/compare/output.xyz"


def test_xyzTranslate_move():
    x: dict = {"file": in_file,
               "move": [5, 0, 0], "out": out_file, "cut": None}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    compare = "tests/compare/xyzTranslate-move.xyz"
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)


def test_xyzTranslate_cut():
    x: dict = {"file": in_file,
               "move": [5, 0, 0], "out": out_file, "cut": 10}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    compare = "tests/compare/xyzTranslate-cut.xyz"
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)


def test_xyzTranslate_cut_move():
    x: dict = {"file": in_file,
               "move": [5, 0, 0], "out": out_file, "cut": 3}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)

    x: dict = {"file": out_file,
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
    x: dict = {"file": in_file, "out": out_file, "cut": 10}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(args)
    assert e.type == SystemExit
    assert e.value.code == 1    # try and exception


def test_xyzTranslate_miss_file():
    x: dict = {"file": "tests/data/crest_conformers000.xyz",
               "out": out_file, "cut": 10}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(args)
    assert e.type == SystemExit
    assert e.value.code == 1    # try and exception


if __name__ == "__main__":
    import cProfile
    cProfile.run("test_xyzTranslate_move()", sort="cumtime")
