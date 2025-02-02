#!/usr/bin/env python3
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzSerial as xyzSerial
import filecmp
import os


def test_xyzSerial_new():
    x: dict = {"file": "tests/data/crest_conformers.xyz", "new": True,
               "keep": False, "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)
    xyzSerial.main(args)

    compare = "tests/compare/xyzSerial-new.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    os.remove(args.out)


def test_xyzSerial_keep():
    x: dict = {"file": "tests/data/crest_conformers_xyzSerial.xyz", "new": False,
               "keep": True, "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)
    xyzSerial.main(args)

    compare = "tests/compare/xyzSerial-keep.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    os.remove(args.out)


def test_xyzSerial_filename_miss():
    x: dict = {"file": "tests/data/crest_conformers_xyzSerial000.xyz", "new": True,
               "keep": True, "out": "tests/compare/output.xyz", "print": False}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        xyzSerial.main(args)
    assert e.type == SystemExit
    assert e.value.code == 1


def test_xyzSerial_miss():
    x: dict = {}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        xyzSerial.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2    # for argprarse wrong
