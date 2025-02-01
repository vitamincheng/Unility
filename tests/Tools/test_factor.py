
import pathlib
import argparse
from pathlib import Path
import pytest
from icecream import ic
from censo_ext.Tools.factor import method_factor_analysis, method_factor_opt
import os
import numpy as np


def test_factor_analysis():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ.xyz", "factor": 0.50}
    args = argparse.Namespace(**x)
    idx1_minor_factor, Table_S = method_factor_analysis(args)
    assert idx1_minor_factor == [np.int64(1), np.int64(
        2), np.int64(7), np.int64(8), np.int64(15), np.int64(19), np.int64(21), np.int64(23), np.int64(24), np.int64(27), np.int64(30), np.int64(33), np.int64(35), np.int64(37), np.int64(40), np.int64(44), np.int64(45), np.int64(50), np.int64(52), np.int64(55)]
    assert len(Table_S) == 29
    assert Table_S[np.int64(1)] == pytest.approx(
        np.float64(0.41868188431231546))
    assert Table_S[np.int64(61)] == pytest.approx(
        np.float64(3.8871420948650774))

    x: dict = {"file": "tests/data/crest_conformers1.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ.xyz", "factor": 0.50}
    args = argparse.Namespace(**x)
    idx1_minor_factor, Table_S = method_factor_analysis(args)
    assert idx1_minor_factor == []
    assert len(Table_S) == 29
    assert Table_S[np.int64(1)] == np.float64(0.0)
    assert Table_S[np.int64(61)] == np.float64(0.0)


def test_method_factor_opt():

    x: dict = {"file": "tests/data/crest_conformers.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ.xyz", "factor": 0.50}
    args = argparse.Namespace(**x)
    idx1_minor_factor, Table_S = method_factor_analysis(args)

    a = method_factor_opt(args, idx1_minor_factor, Table_S)
    assert a[0] == True
    assert a[1] == [np.int64(52), 55]
    assert a[2] == pytest.approx(np.float64(0.9863091021443952))
