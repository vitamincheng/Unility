#!/usr/bin/env python3
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzSerial as xyzSerial


def test_xyzSerial_new():
    x: dict = {"file": "tests/data/crest_conformers.xyz", "new": True,
               "keep": False, "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)
    xyzSerial.main(args)

    import filecmp
    compare = "tests/compare/xyzSerial-new.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    import os
    os.remove(args.out)


def test_xyzSerial_keep():
    x: dict = {"file": "tests/data/crest_conformers_xyzSerial.xyz", "new": False,
               "keep": True, "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)
    xyzSerial.main(args)

    import filecmp
    compare = "tests/compare/xyzSerial-keep.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    import os
    os.remove(args.out)
