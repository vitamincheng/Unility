import argparse
from pathlib import Path
import pytest
import numpy as np
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.Tools.calculate_rmsd import cal_RMSD_xyz
from censo_ext.Tools.calculate_rmsd import atom2str, atom2int, rmsd, kabsch_rotate, kabsch, centroid, get_Coordinates


def test_calculate_rmsd_miss_args():

    xyzFile: GeometryXYZs = GeometryXYZs()

    with pytest.raises(FileNotFoundError) as e:
        xyzFile.method_read_xyz()
    assert str(e.value) == ". was not found or is a directory"


def test_calculate_rmsd_all():
    # for crest_conformers.xyz
    xyzFile: GeometryXYZs = GeometryXYZs(
        Path("tests/data/crest_conformers.xyz"))
    xyzFile.method_read_xyz()
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": False, "debug": False}
    idx_atom1 = cal_RMSD_xyz(xyzFile, 1, 2, args=argparse.Namespace(**x))
    assert len(idx_atom1[0]) == 73
    assert idx_atom1[1] == pytest.approx(1.7608313888081009)

    x = {"remove_idx": None, "add_idx": None, "bond_broken": [
        52, 55], "ignore_Hydrogen": True, "debug": False}
    idx_atom1 = cal_RMSD_xyz(xyzFile, 1, 2, args=argparse.Namespace(**x))
    assert len(idx_atom1[0]) == 24
    assert idx_atom1[0][1] == pytest.approx(0.0022804570676915915)
    assert idx_atom1[0][24] == pytest.approx(0.02272775019255817)
    assert idx_atom1[1] == pytest.approx(0.35959367835047723)

    # for isomers.xyz
    xyzFile = GeometryXYZs(
        Path("tests/data/isomers.xyz"))
    xyzFile.method_read_xyz()
    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": False, "debug": False}
    idx_atom1 = cal_RMSD_xyz(xyzFile, 1, 2, args=argparse.Namespace(**x))

    assert len(idx_atom1[0]) == 17
    assert idx_atom1[0][1] == pytest.approx(0.052818544187445145)
    assert idx_atom1[0][17] == pytest.approx(0.5847626056423644)
    assert idx_atom1[1] == pytest.approx(0.9657626138106812)


def test_calculate_atom2str():
    assert atom2str(1) == "H"
    assert atom2str(8) == "O"
    assert atom2str(6) == "C"


def test_calculate_atom2int():
    assert atom2int("C") == 6
    assert atom2int("O") == 8
    assert atom2int("N") == 7


def test_calculate_rmsd():
    P = np.array([[0, 0], [3, 4]])
    Q = np.array([[2, 3], [4, 5]])
    idx_atom = [0, 1]
    expected_square: dict[int, float] = {0: 13.0, 1: 2.0}
    expected_rmsd: float = 2.7386127875
    square, rmsd_value = rmsd(P, Q, idx_atom)
    assert square == expected_square
    assert np.isclose(rmsd_value, expected_rmsd)


def test_calculate_kabsch_rotate():
    P = np.array([[0, 0], [3, 4]])
    Q = np.array([[2, 3], [4, 5]])
    U = kabsch(P, Q)
    Res = kabsch_rotate(P, Q)
    assert np.allclose(Res, np.dot(P, U))


def test_calculate_kabsch():
    P = np.array([[0, 0], [3, 4]])
    Q = np.array([[2, 3], [4, 5]])
    expected_U = np.array([
        [0.9995120, -0.0312347],
        [0.0312347,  0.9995120]
    ])
    U = kabsch(P, Q)
    assert np.allclose(U, expected_U)


def test_calculate_centroid():
    V = np.array([[0, 0, 0], [0, 0, 1], [1, 0, 0]])
    expected_c = np.array([0.33333333, 0.00, 0.33333333])
    c = centroid(V)
    assert np.allclose(c, expected_c)


def test_calculate_get_Coordinates():
    xyzFile = GeometryXYZs("tests/data/isomers.xyz")
    xyzFile.method_read_xyz()
    element, V = get_Coordinates(xyzFile, 0)
    assert len(element) == 17
    assert V[0][0] == pytest.approx(1.8513508441)
    assert V[0][2] == pytest.approx(-0.0716885123)
    # assert element == expected_element
    # assert np.allclose(V, expected_V)


def test_calculate_cal_RMSD_xyz():
    xyzFile = GeometryXYZs("tests/data/isomers.xyz")
    xyzFile.method_read_xyz()
    idx_p = 1
    idx_q = 2
    # expected_square: dict[int, float] = {0: 5.0}
    # expected_rmsd: float = np.sqrt(5.0 / 3)

    x: dict = {"remove_idx": None, "add_idx": None, "bond_broken": None,
               "ignore_Hydrogen": False, "debug": False}
    square, rmsd_value = cal_RMSD_xyz(
        xyzFile, idx_p, idx_q, args=argparse.Namespace(**x))
    assert len(square) == 17
    assert np.isclose(rmsd_value, 0.9657626138106813)
