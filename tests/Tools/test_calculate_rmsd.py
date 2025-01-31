import pathlib
import argparse
from pathlib import Path
import pytest
import numpy as np
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.Tools.calculate_rmsd import cal_rmsd_xyz


def test_calculate_rmsd():
    # for crest_conformers.xyz
    xyzfile = GeometryXYZs(
        Path("tests/data/crest_conformers.xyz"))
    xyzfile.method_read_xyz()
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": False, "debug": False}
    idx_atom1 = cal_rmsd_xyz(xyzfile, 1, 2, args=argparse.Namespace(**x))
    assert len(idx_atom1[0]) == 73
    assert idx_atom1[1] == np.float64(1.7608313888081009)

    x = {"remove_idx": None, "add_idx": None, "bond_broken": [
        52, 55], "ignore_Hydrogen": True, "debug": False}
    idx_atom1 = cal_rmsd_xyz(xyzfile, 1, 2, args=argparse.Namespace(**x))
    assert len(idx_atom1[0]) == 24
    assert idx_atom1[1] == pytest.approx(np.float64(0.35959367835047723))

    # for isomers.xyz
    xyzfile = GeometryXYZs(
        Path("tests/data/isomers.xyz"))
    xyzfile.method_read_xyz()
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": False, "debug": False}
    idx_atom1 = cal_rmsd_xyz(xyzfile, 1, 2, args=argparse.Namespace(**x))

    assert len(idx_atom1[0]) == 17
    assert idx_atom1[1] == pytest.approx(np.float64(0.9657626138106812))
