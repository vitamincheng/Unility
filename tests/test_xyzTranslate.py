#!/usr/bin/env python
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzTranslate as xyzTranslate
import filecmp
import os
inFile: Path = Path("tests/data/crest_conformers.xyz")
outFile: Path = Path("tests/compare/output.xyz")


def test_xyzTranslate_move():
    x: dict = {"file": inFile,
               "move": [5, 0, 0], "out": outFile, "cut": None}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    compare: Path = Path("tests/compare/xyzTranslate-move.xyz")
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


def test_xyzTranslate_cut():
    x: dict = {"file": inFile,
               "move": [5, 0, 0], "out": outFile, "cut": 10}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)
    compare: Path = Path("tests/compare/xyzTranslate-cut.xyz")
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


def test_xyzTranslate_cut_move():
    x: dict = {"file": inFile,
               "move": [5, 0, 0], "out": outFile, "cut": 3}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)

    x: dict = {"file": outFile,
               "move": [0, 0, 5], "out": Path("tests/compare/output2.xyz"), "cut": 3}
    args = argparse.Namespace(**x)
    xyzTranslate.main(args)

    compare: Path = Path("tests/compare/xyzTranslate-cut-move.xyz")
    assert filecmp.cmp(args.out, compare)
    os.remove(args.file)
    os.remove(args.out)


def test_xyzTranslate_miss_args():
    x: dict = {}
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_xyzTranslate_without_move():
    x: dict = {"file": inFile, "out": outFile, "cut": 10}
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 1    # try and exception


def test_xyzTranslate_miss_file():
    x: dict = {"file": Path("tests/data/crest_conformers000.xyz"),
               "out": outFile, "cut": 10}
    with pytest.raises(SystemExit) as e:
        xyzTranslate.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 1    # try and exception


# if __name__ == "__main__":
#    import cProfile
#    cProfile.run("test_xyzTranslate_move()", sort="cumtime")
