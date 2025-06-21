#!/usr/bin/env python3
import os
import pytest
from pathlib import Path
import argparse
import censo_ext.molclus_thermo as thermo
import filecmp
import platform

in_file = f"tests/data/EthylAcetate/01.Crest/crest_conformers.xyz"


def test_molclus_thermo_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        thermo.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_molclus_thermo_alpb():
    x: dict = {"file": in_file, "alpb": "CHCl3", "gbsa": None, "chrg": 0,
               "uhf": 1, "opt": True, "method": "gfn2"}
    args = argparse.Namespace(**x)

    Result = thermo.main(args)

    assert float(Result[0]) == pytest.approx(0.083147421, abs=0.0000002)
    assert float(Result[1]) == pytest.approx(0.082439908, abs=0.0000002)

    import subprocess
    subprocess.call(
        f"rm -f charges g98.out hessian thermo.out vibspectrum xtbhess.xyz", shell=True)
    subprocess.call(
        f"rm -f wbo xtb_enso.json xtbopt.xyz xtbrestart xtbtopo.mol", shell=True)
