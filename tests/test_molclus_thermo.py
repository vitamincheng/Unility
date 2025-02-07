#!/usr/bin/env python3
import os
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_thermo as thermo
import filecmp
import platform


def test_molclus_thermo_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        thermo.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_molclus_thermo_alpb():
    x: dict = {"file": "tests/data/crest_conformers2.xyz", "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1, "opt": True,
               "method": "gfn2"}
    args = argparse.Namespace(**x)

    Result = thermo.main(args)

    assert float(Result[0]) == pytest.approx(0.5922643, abs=0.0000002)
    assert float(Result[1]) == pytest.approx(0.5926563, abs=0.0000002)

    import subprocess
    subprocess.call(
        f"rm -f charges g98.out hessian thermo.out vibspectrum", shell=True)
    subprocess.call(
        f"rm -f wbo xtb_enso.json xtbopt.xyz xtbrestart xtbtopo.mol", shell=True)
