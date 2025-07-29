import argparse
from pathlib import Path
import pytest
from censo_ext.Tools.topo import Topo


def test_topo_FileName_miss_args():
    x = {"file": Path("test.xyz"),
         "bonding": 20, "print": True, "debug": False}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        Sts_topo: Topo = Topo(Path(args.file))
    assert e.type == SystemExit
    assert e.value.code == 1  # for argparse error


@pytest.mark.parametrize(argnames="input_Path,bonding,CN_Dict",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"), 51,
                                     {1: 3, 2: 4, 3: 4, 4: 1, 5: 1, 6: 3, 7: 4, 8: 4, 9: 1, 10: 1, 11: 1, 12: 1, 13: 2, 14: 1, 15: 3, 16: 3, 17: 1, 18: 1, 19: 3, 20: 1, 21: 3, 22: 1, 23: 4, 24: 4, 25: 1, 26: 1, 27: 4, 28: 1, 29: 1, 30: 4, 31: 1, 32: 1, 33: 4, 34: 1, 35: 4, 36: 1,
                                      37: 4, 38: 1, 39: 1, 40: 4, 41: 1, 42: 1, 43: 1, 44: 4, 45: 4, 46: 1, 47: 1, 48: 1, 49: 1, 50: 4, 51: 3, 52: 3, 53: 1, 54: 1, 55: 4, 56: 4, 57: 4, 58: 1, 59: 1, 60: 4, 61: 4, 62: 1, 63: 1, 64: 1, 65: 1, 66: 1, 67: 1, 68: 1, 69: 1, 70: 1, 71: 1, 72: 1, 73: 1}),
                                    (Path("tests/data/isomers.xyz"), 3,
                                     {1: 4, 2: 4, 3: 4, 4: 1, 5: 1, 6: 3, 7: 1, 8: 1,
                                        9: 4, 10: 4, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1})])
def test_topo_get_cn(input_Path: str, bonding: int, CN_Dict: dict):
    x = {"file": input_Path,
         "bonding": bonding, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    assert Sts_topo.get_cn() == CN_Dict


@pytest.mark.parametrize(argnames="input_Path,bonding,Neighbors_Atoms",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"), 51, [44, 52]),
                                    (Path("tests/data/isomers.xyz"), 3, [2, 6])])
def test_topo_Bonding(input_Path: Path, bonding: int, Neighbors_Atoms: list):
    # for crest_conformers.xyz
    x = {"file": input_Path,
         "bonding": bonding, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    assert Sts_topo.method_bonding(args) == Neighbors_Atoms


@pytest.mark.parametrize(argnames="input_Path,bonding,len_neighbor,circle_Mols,residual_Mols",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"), 51, 29,
                                     [[1, 6, 7, 8, 3, 2], [21, 33, 30, 27, 24, 23], [
                                         21, 23, 24, 27, 30, 40, 37, 35, 33], [33, 30, 40, 37, 35]],
                                     [[3, 13], [1, 15, 19, 21], [6, 16], [
                                         40, 44, 50, 51, 52, 55, 56, 57, 60, 61], [30, 45]]
                                     ),
                                    (Path("tests/data/isomers.xyz"), 3, 7,
                                     [[6, 3, 2, 1, 10, 9]], [[6, 17]],
                                     )])
def test_topo_topology(input_Path: Path, bonding: int, len_neighbor: int, circle_Mols: list, residual_Mols: list):
    # for crest_conformers.xyz
    x = {"file": input_Path, "bonding": bonding, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    mol, neighbors, circle_Mols_R, residual_Mols_R = Sts_topo.topology()
    assert len(neighbors) == len_neighbor
    assert (circle_Mols_R) == circle_Mols
    assert residual_Mols_R == residual_Mols


@pytest.mark.parametrize(argnames="input_Path,bond_broken,broken_Atoms_H",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"), [40, 44],
                                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                                      24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 48]),
                                    (Path("tests/data/isomers.xyz"), [3, 2],
                                     [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])])
def test_topo_Broken_bond_H(input_Path: Path, bond_broken: list, broken_Atoms_H: list):
    # for crest_conformers.xyz
    x = {"file": input_Path, "bond_broken": bond_broken,
         "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    assert Sts_topo.method_broken_bond_H(args) == broken_Atoms_H


@pytest.mark.parametrize(argnames="input_Path,bond_broken,broken_Atoms",
                         argvalues=[(Path("tests/data/crest_conformers.xyz"), [52, 55],
                                     [1, 2, 3, 6, 7, 8, 13, 15, 16, 19, 21,
                                      23, 24, 27, 30, 33, 35, 37, 40, 44, 45, 50, 51, 52]),
                                    (Path("tests/data/isomers.xyz"), [3, 2],
                                        [1, 3, 6, 9, 10, 17])])
def test_topo_Broken_bond(input_Path: Path, bond_broken: list, broken_Atoms: list):
    x = {"file": input_Path,
         "bond_broken": bond_broken, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    assert Sts_topo.method_broken_bond(args) == broken_Atoms
