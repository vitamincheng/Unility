#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.xyzGenFlexible as xyzGenFlexible
import filecmp
import platform


def test_xyzGenflexible_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        xyzGenFlexible.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_xyzGenFlexible_args():
    x: dict = {"file": "tests/data/crest_conformers.xyz", "manual": False,
               "out": "tests/compare/xyzGen_isomers.xyz"}
    args = argparse.Namespace(**x)
    xyzGenFlexible.main(args)

    compare = "tests/compare/xyzGen_Darwin.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    os.remove(args.out)


def test_xyzGenFlexible_args_manual(monkeypatch):
    x: dict = {"file": "tests/data/crest_conformers.xyz", "manual": True,
               "out": "tests/compare/xyzGen_isomers.xyz"}
    args = argparse.Namespace(**x)
    import io
    monkeypatch.setattr('sys.stdin', io.StringIO("55"))
    xyzGenFlexible.main(args)

    compare = "tests/compare/xyzGen_Darwin_manual.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)
    os.remove(args.out)
