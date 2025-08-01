#!/usr/bin/env python3
import pytest
import argparse
import censo_ext.molclus_thermo as thermo

inFile = "tests/data/06.EthylAcetate/01.Crest/crest_conformers.xyz"


def test_molclus_thermo_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        thermo.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_molclus_thermo_alpb():
    x: dict = {"file": inFile, "alpb": "CHCl3", "gbsa": None, "chrg": 0,
               "uhf": 1, "opt": True, "method": "gfn2"}
    args = argparse.Namespace(**x)

    Res: list[str] = thermo.main(args)

    assert float(Res[0]) == pytest.approx(0.083147421, abs=0.0000002)
    assert float(Res[1]) == pytest.approx(0.082439908, abs=0.0000002)

    import subprocess
    subprocess.call(
        "rm -f charges g98.out hessian thermo.out vibspectrum xtbhess.xyz", shell=True)
    subprocess.call(
        "rm -f wbo xtb_enso.json xtbopt.xyz xtbrestart xtbtopo.mol", shell=True)
