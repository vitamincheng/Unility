#!/usr/bin/env python
import pytest
import argparse
import censo_ext.xyzRotate as xyzRotate
import os
import filecmp
import platform
from pathlib import Path
_system = platform.system()
inFile: Path = Path("tests/data/crest_conformers.xyz")


def test_xyzRotate_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzRotate.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_xyzRotate_3_1():
    x: dict = {"file": inFile, "atoms": [52, 55], "cuts": 3, "spec": 1,
               "print": False, "check": False,
               "out": Path("tests/compare/output_xyzReturnOandZ.xyz")}
    args = argparse.Namespace(**x)
    xyzRotate.main(args)

    compare: Path = Path("tests/compare/xyzRotate_3_1.xyz")
    assert filecmp.cmp(args.out, compare)  # type: ignore
    os.remove(args.out)


def test_xyzRotate_3_2():
    x: dict = {"file": inFile, "atoms": [52, 55], "cuts": 3, "spec": 2,
               "print": False, "check": False,
               "out": Path("tests/compare/output_xyzReturnOandZ.xyz")}
    args = argparse.Namespace(**x)
    xyzRotate.main(args)

    compare: Path = Path("tests/compare/xyzRotate_3_2.xyz")
    assert filecmp.cmp(args.out, compare)  # type: ignore
    os.remove(args.out)
