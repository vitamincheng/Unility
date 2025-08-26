#!/usr/bin/env python
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzSerial as xyzSerial
import filecmp
import os
outFile: Path = Path("tests/compare/output.xyz")


def test_xyzSerial_new():
    x: dict = {"file": "tests/data/crest_conformers.xyz", "new": True,
               "keep": False, "out": outFile, "print": False}
    args = argparse.Namespace(**x)
    xyzSerial.main(args)

    compare = "tests/compare/xyzSerial-new.xyz"
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


def test_xyzSerial_keep():
    x: dict = {"file": "tests/data/crest_conformers4.xyz", "new": False,
               "keep": True, "out": outFile, "print": False}
    args = argparse.Namespace(**x)
    xyzSerial.main(args)

    compare = "tests/compare/xyzSerial-keep.xyz"
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


def test_xyzSerial_filename_miss():
    x: dict = {"file": "tests/data/crest_conformers_xyzSerial000.xyz", "new": True,
               "keep": True, "out": outFile, "print": False}

    args = argparse.Namespace(**x)

    with pytest.raises(FileNotFoundError) as e:
        xyzSerial.main(args)
    assert str(e.value) == f"{x['file']} The file is not Exist ..."


def test_xyzSerial_miss():
    x: dict = {}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        xyzSerial.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2    # for argprarse wrong
