#!/usr/bin/env python
import pytest
import argparse
# import filecmp
import censo_ext.ensoGenFlexible as ensoGenFlexible
from censo_ext.Tools.utility import delete_all_files


def test_ensoGenFlexible_miss_args():
    x: dict = {}
    with pytest.raises(SystemExit) as e:
        ensoGenFlexible.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


@pytest.mark.slow
def test_ensoGenFlexible_single_xyz() -> None:
    inFile = "tests/data/06.EthylAcetate/traj.xyz"
    x: dict = {"file": inFile, "manual": False, "temp": 298.15}
    ensoGenFlexible.main(argparse.Namespace(**x))
    # compare = Path("tests/compare/molManipulate.xyz")
    # assert filecmp.cmp(args.out, compare)
    delete_all_files("anmr_enso.new")


@pytest.mark.slow
def test_ensoGenFlexible_multi_xyzs() -> None:
    inFile = "tests/data/06.EthylAcetate/03.Censo/crest_conformers.xyz"
    x: dict = {"file": inFile, "manual": False, "temp": 298.15}
    ensoGenFlexible.main(argparse.Namespace(**x))
    # compare = Path("tests/compare/molManipulate.xyz")
    # assert filecmp.cmp(args.out, compare)
    delete_all_files("anmr_enso.new")
