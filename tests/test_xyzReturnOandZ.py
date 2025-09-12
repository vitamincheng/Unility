#!/usr/bin/env python
import pytest
import argparse
import censo_ext.xyzReturnOandZ as xyzReturnOandZ
import os
import filecmp
import platform
from pathlib import Path
_system = platform.system()
inFile: Path = Path("tests/data/crest_conformers.xyz")


def test_xyzReturnOandZ_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzReturnOandZ.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_xyzReturnOandZ():
    x: dict = {"file": inFile,
               "atom": [30, 45, 47], "print": False, "replace": False, "out": Path("tests/compare/output_xyzReturnOandZ.xyz")}
    args = argparse.Namespace(**x)
    xyzReturnOandZ.main(args)

    if _system == 'Darwin':
        compare: Path = Path("tests/compare/xyzReturnOandZ_Darwin.xyz")
    elif _system == "Linux":
        compare: Path = Path("tests/compare/xyzReturnOandZ.xyz")
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare)  # type: ignore
    os.remove(args.out)


def test_xyzReturnOandZ_auto():
    x: dict = {"file": inFile, "auto": True,
               "atom": None, "print": False, "replace": False, "out": Path("tests/compare/output_xyzReturnOandZ_auto.xyz")}
    args = argparse.Namespace(**x)
    xyzReturnOandZ.main(args)

    if _system == 'Darwin':
        compare: Path = Path("tests/compare/xyzReturnOandZ_auto_Darwin.xyz")
    elif _system == "Linux":
        compare: Path = Path("tests/compare/xyzReturnOandZ_auto.xyz")
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare)  # type: ignore
    os.remove(args.out)
