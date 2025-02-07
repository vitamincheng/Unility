#!/usr/bin/env python3
import os
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_orca as orca
import filecmp
import platform


def test_orca_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        orca.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_orca_opt():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "template": "template.inp",
               "remove": True, "out": "tests/compare/orca_isomers.xyz"}

    args = argparse.Namespace(**x)
    if platform.system() != "Darwin":
        orca.main(args)
