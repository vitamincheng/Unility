#!/usr/bin/env python
import pytest
from pathlib import Path
import argparse
import censo_ext.FactorAnalysis as FactorAnalysis
inFile: Path = Path("tests/data/crest_conformers.xyz")


def test_FactorAnalysis_miss_args():
    x: dict = {}
    with pytest.raises(SystemExit) as e:
        FactorAnalysis.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_FactorAnalysis_opt():
    x: dict = {"file": inFile, "factor": 0.50, "opt": True, "thr": 2,
               "ignore_Hydrogen": False, "Analysis": True, "Filter": None}
    assert FactorAnalysis.main(argparse.Namespace(**x)) is None


def test_FactorAnalysis_A_bond_broken() -> None:
    x: dict = {"file": inFile, "factor": 0.50, "opt": False, "thr": 2,
               "bond_broken": [40, 44], "ignore_Hydrogen": False,
               "Analysis": True, "Filter": None}
    assert FactorAnalysis.main(argparse.Namespace(**x)) is None


def test_FactorAnalysis_F_bond_broken() -> None:
    x: dict = {"file": inFile, "factor": 0.50, "opt": False, "thr": 2,
               "bond_broken": [40, 44], "ignore_Hydrogen": True, "Analysis": None,
               "remove_idx": None, "add_idx": None, "Filter": True}
    out_print = "result.txt"
    import sys
    o_stdout = sys.stdout
    with open(out_print, "w") as f:
        sys.stdout = f
        FactorAnalysis.main(argparse.Namespace(**x))
    sys.stdout = o_stdout
    with open(out_print, "r") as f:
        lines: list[str] = f.readlines()

    assert float(lines[-1].split()[-1]) == pytest.approx(1.0651832633943255)

    import subprocess
    subprocess.call("rm -rf Final_Result", shell=True)
    subprocess.call(f"rm -rf {out_print}", shell=True)


def test_FactorAnalysis_F_bond_broken_incl_H() -> None:
    x: dict = {"file": inFile, "factor": 0.50, "opt": False, "thr": 2,
               "bond_broken": [40, 44], "ignore_Hydrogen": False, "Analysis": None,
               "remove_idx": None, "add_idx": None, "Filter": True}

    with pytest.raises(ValueError) as e:
        FactorAnalysis.main(argparse.Namespace(**x))
    assert str(e.value) == "Only support under ignore Hydrogen condition "
