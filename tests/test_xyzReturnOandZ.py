#!/usr/bin/env python3
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzReturnOandZ as xyzReturnOandZ
import os
import filecmp
import platform
_system = platform.system()
inFile = f"tests/data/crest_conformers.xyz"


def test_xyzReturnOandZ_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzReturnOandZ.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_xyzReturnOandZ():
    x: dict = {"file": inFile,
               "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ.xyz"}
    args = argparse.Namespace(**x)
    xyzReturnOandZ.main(args)

    compare = f""
    if _system == 'Darwin':
        compare = f"tests/compare/xyzReturnOandZ_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/xyzReturnOandZ.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)


def test_xyzReturnOandZ_auto():
    x: dict = {"file": inFile, "auto": True,
               "atom": None, "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ_auto.xyz"}
    args = argparse.Namespace(**x)
    xyzReturnOandZ.main(args)

    compare = ""
    if _system == 'Darwin':
        compare = f"tests/compare/xyzReturnOandZ_auto_Darwin.xyz"
    elif _system == "Linux":
        compare = f"tests/compare/xyzReturnOandZ_auto.xyz"
    else:
        pytest.raises(
            ValueError, match="OS system only can run under Darwin or Linux")
    assert filecmp.cmp(args.out, compare) == True
    os.remove(args.out)
