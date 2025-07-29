#!/usr/bin/env python3
import pytest
from pathlib import Path
import argparse
import censo_ext.datNormalized as datNormalized
import filecmp
from censo_ext.Tools.utility import delete_all_files
import platform


def test_datNormalized_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        datNormalized.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_datNormalized():
    in_file: Path = Path("tests/data/04.Hydrogen/anmr.dat")
    out_file: Path = Path("tests/compare/out.dat")
    compare_file: Path = Path("tests/compare/anmr_normal.dat")

    x: dict = {"file": in_file, "start": -5.0, "end": 15.0, "dpi": 10000,
               "out": out_file}
    args = argparse.Namespace(**x)
    datNormalized.main(args)
    assert filecmp.cmp(out_file, compare_file) == True
    delete_all_files(out_file)
