#!/usr/bin/env python
import filecmp
import pytest
import argparse
import censo_ext.molclus_thermo as thermo
from pathlib import Path

from censo_ext.Tools.utility import delete_all_files

inFile: Path = Path("tests/data/06.EthylAcetate/01.Crest/crest_conformers.xyz")


def test_molclus_thermo_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        thermo.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_molclus_thermo_alpb():
    x: dict = {"file": inFile, "alpb": "CHCl3", "gbsa": None, "chrg": 0,
               "uhf": 1, "opt": True, "method": "gfn2", "enso": True, "temp": 298.15}

    thermo.main(argparse.Namespace(**x))
    target = Path("anmr_enso.new")
    compare = Path("tests/compare/thermo/anmr_enso.new")
    assert filecmp.cmp(compare, target)
    delete_all_files(target)
