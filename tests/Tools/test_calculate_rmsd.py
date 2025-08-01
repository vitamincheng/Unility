import argparse
from pathlib import Path
import pytest
# import numpy as np
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz


def test_calculate_rmsd_miss_args():

    xyzfile = GeometryXYZs()

    with pytest.raises(FileNotFoundError) as e:
        xyzfile.method_read_xyz()
    assert str(e.value) == ". was not found or is a directory"


def test_calculate_rmsd():
    # for crest_conformers.xyz
    xyzfile = GeometryXYZs(
        Path("tests/data/crest_conformers.xyz"))
    xyzfile.method_read_xyz()
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": False, "debug": False}
    idx_atom1 = cal_RMSD_xyz(xyzfile, 1, 2, args=argparse.Namespace(**x))
    assert len(idx_atom1[0]) == 73
    assert idx_atom1[1] == pytest.approx(1.7608313888081009)

    x = {"remove_idx": None, "add_idx": None, "bond_broken": [
        52, 55], "ignore_Hydrogen": True, "debug": False}
    idx_atom1 = cal_RMSD_xyz(xyzfile, 1, 2, args=argparse.Namespace(**x))
    assert len(idx_atom1[0]) == 24
    assert idx_atom1[0][1] == pytest.approx(0.0022804570676915915)
    assert idx_atom1[0][24] == pytest.approx(0.02272775019255817)
    assert idx_atom1[1] == pytest.approx(0.35959367835047723)

    # for isomers.xyz
    xyzfile = GeometryXYZs(
        Path("tests/data/isomers.xyz"))
    xyzfile.method_read_xyz()
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": False, "debug": False}
    idx_atom1 = cal_RMSD_xyz(xyzfile, 1, 2, args=argparse.Namespace(**x))

    assert len(idx_atom1[0]) == 17
    assert idx_atom1[0][1] == pytest.approx(0.052818544187445145)
    assert idx_atom1[0][17] == pytest.approx(0.5847626056423644)
    assert idx_atom1[1] == pytest.approx(0.9657626138106812)
