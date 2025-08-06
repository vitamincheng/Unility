#!/usr/bin/env python
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzGenFlexible as xyzGenFlexible
import filecmp

inFile: Path = Path("tests/data/crest_conformers.xyz")
outFile: Path = Path("tests/compare/xyzGen_isomers.xyz")


def test_xyzGenflexible_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzGenFlexible.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_xyzGenFlexible_args():
    x: dict = {"file": inFile, "manual": False,
               "out": outFile}
    args = argparse.Namespace(**x)
    xyzGenFlexible.main(args)

    compare = "tests/compare/xyzGen_Darwin.xyz"
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)


def test_xyzGenFlexible_args_manual(monkeypatch):
    x: dict = {"file": inFile, "manual": True,
               "out": outFile}
    args = argparse.Namespace(**x)
    import io
    monkeypatch.setattr('sys.stdin', io.StringIO("55"))
    xyzGenFlexible.main(args)

    compare = "tests/compare/xyzGen_Darwin_manual.xyz"
    assert filecmp.cmp(args.out, compare)
    os.remove(args.out)
