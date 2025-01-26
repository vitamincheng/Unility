#!/usr/bin/env python3
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzReturnOandZ as xyzReturnOandZ


def test_xyzReturnOandZ():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output.xyz"}
    args = argparse.Namespace(**x)
    xyzReturnOandZ.main(args)

    import filecmp
    import platform
    print(platform.system())
    if platform.system() == 'Darwin':
        compare = "tests/compare/xyzReturnOandZ_Darwin.xyz"
    else:
        compare = "tests/compare/xyzReturnOandZ.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    import os
    os.remove(args.out)


def test_xyzReturnOandZ_auto():
    x: dict = {"file": "tests/data/crest_conformers.xyz", "auto": True,
                       "atom": None, "print": False, "replace": False, "out": "tests/compare/output1.xyz"}
    args = argparse.Namespace(**x)
    xyzReturnOandZ.main(args)

    import filecmp
    import platform
    if platform.system() == 'Darwin':
        compare = "tests/compare/xyzReturnOandZ-auto_Darwin.xyz"
    else:
        compare = "tests/compare/xyzReturnOandZ-auto.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    # import os
    # os.remove(args.out)
