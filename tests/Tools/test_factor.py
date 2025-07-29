
import argparse
from pathlib import Path
import pytest
from icecream import ic
from censo_ext.Tools.factor import method_factor_analysis, method_factor_opt
import numpy as np


def test_factor_analysis_miss_args():
    x: dict = {"file": "tests/data/crest_conformers000.xyz", "print": False,
               "replace": False, "factor": 0.5}
    args = argparse.Namespace(**x)

    with pytest.raises(FileNotFoundError) as e:
        idx1_minor_factor, Table_S = method_factor_analysis(args)
    # for argparse error
    assert str(e.value) == f"{x["file"]} The file is not Exist ..."


def test_factor_analysis():
    x: dict = {"file": "tests/data/crest_conformers.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ.xyz", "factor": 0.50}
    args = argparse.Namespace(**x)
    idx1_minor_factor, Table_S = method_factor_analysis(args)
    assert idx1_minor_factor == [1, 2, 7, 8, 15, 19, 21,
                                 23, 24, 27, 30, 33, 35, 37, 40, 44, 45, 50, 52, 55]
    assert len(Table_S) == 29
    assert Table_S[1] == pytest.approx(0.41868188431231546)
    assert Table_S[61] == pytest.approx(3.8871420948650774)

    x: dict = {"file": "tests/data/crest_conformers1.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ.xyz", "factor": 0.50}
    args = argparse.Namespace(**x)
    idx1_minor_factor, Table_S = method_factor_analysis(args)
    assert idx1_minor_factor == []
    assert len(Table_S) == 29
    assert Table_S[1] == 0.0
    assert Table_S[61] == 0.0


def test_factor_opt():

    x: dict = {"file": "tests/data/crest_conformers.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "tests/compare/output_xyzReturnOandZ.xyz", "factor": 0.50}
    args = argparse.Namespace(**x)
    idx1_minor_factor, Table_S = method_factor_analysis(args)

    a0, a1, a2 = method_factor_opt(
        args, idx1_minor_factor, Table_S)  # type: ignore
    assert a0 == True
    assert a1 == [52, 55]
    assert a2 == pytest.approx(0.9863091021443952)
