#!/usr/bin/env python3
import os
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.datNormalized as datNormalized
import filecmp
import platform


def test_datNormalized_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        datNormalized.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_datNormalized():
    in_file = "tests/data/04.Hydrogen/anmr.dat"
    out_file = "tests/compare/out.dat"
    x: dict = {"file": in_file, "start": -5.0, "end": 15.0, "dpi": 10000,
               "out": out_file}
    args = argparse.Namespace(**x)
    datNormalized.main(args)
    assert filecmp.cmp(out_file, "tests/compare/anmr_normal.dat") == True
    os.remove(out_file)
