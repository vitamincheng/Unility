import argparse
from pathlib import Path
import pytest
from censo_ext.Tools.topo import Topo


def test_method_get_cn():
    # for crest_confomrers.xyz
    x = {"file": "tests/data/crest_conformers.xyz",
         "bonding": 51, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    a = Sts_topo.get_cn()
    assert a == {1: 3, 2: 4, 3: 4, 4: 1, 5: 1, 6: 3, 7: 4, 8: 4, 9: 1, 10: 1, 11: 1, 12: 1, 13: 2, 14: 1, 15: 3, 16: 3, 17: 1, 18: 1, 19: 3, 20: 1, 21: 3, 22: 1, 23: 4, 24: 4, 25: 1, 26: 1, 27: 4, 28: 1, 29: 1, 30: 4, 31: 1, 32: 1, 33: 4, 34: 1, 35: 4, 36: 1,
                 37: 4, 38: 1, 39: 1, 40: 4, 41: 1, 42: 1, 43: 1, 44: 4, 45: 4, 46: 1, 47: 1, 48: 1, 49: 1, 50: 4, 51: 3, 52: 3, 53: 1, 54: 1, 55: 4, 56: 4, 57: 4, 58: 1, 59: 1, 60: 4, 61: 4, 62: 1, 63: 1, 64: 1, 65: 1, 66: 1, 67: 1, 68: 1, 69: 1, 70: 1, 71: 1, 72: 1, 73: 1}
    # for isomers.xyz
    x = {"file": "tests/data/isomers.xyz",
         "bonding": 3, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    a = Sts_topo.get_cn()
    assert a == {1: 4, 2: 4, 3: 4, 4: 1, 5: 1, 6: 3, 7: 1, 8: 1,
                 9: 4, 10: 4, 11: 1, 12: 1, 13: 1, 14: 1, 15: 1, 16: 1, 17: 1}


def test_topo_Bonding():
    # for crest_conformers.xyz
    x = {"file": "tests/data/crest_conformers.xyz",
         "bonding": 51, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    Neighbors_Atoms = Sts_topo.method_bonding(args)
    assert Neighbors_Atoms == [44, 52]
    # for isomers.xyz
    x = {"file": "tests/data/isomers.xyz",
         "bonding": 3, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    Neighbors_Atoms = Sts_topo.method_bonding(args)
    assert Neighbors_Atoms == [2, 6]


def test_topo_topology():
    # for crest_conformers.xyz
    x = {"file": "tests/data/crest_conformers.xyz",
         "bonding": 51, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    mol, neighbors, circle_Mols, residual_Mols = Sts_topo.topology()
    assert len(neighbors) == 29
    assert (circle_Mols) == [[1, 6, 7, 8, 3, 2], [21, 33, 30, 27, 24, 23], [
        21, 23, 24, 27, 30, 40, 37, 35, 33], [33, 30, 40, 37, 35]]
    assert residual_Mols == [[3, 13], [1, 15, 19, 21], [6, 16], [
        40, 44, 50, 51, 52, 55, 56, 57, 60, 61], [30, 45]]
    # for isomers.xyz
    x = {"file": "tests/data/isomers.xyz",
         "bonding": 3, "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    mol, neighbors, circle_Mols, residual_Mols = Sts_topo.topology()
    assert len(neighbors) == 7
    assert (circle_Mols) == [[6, 3, 2, 1, 10, 9]]
    assert residual_Mols == [[6, 17]]


def test_topo_Broken_bond_H():
    # for crest_conformers.xyz
    x = {"file": "tests/data/crest_conformers.xyz",
         "bond_broken": [40, 44], "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    Result = Sts_topo.method_broken_bond_H(args)
    assert Result == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                      24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 48]
    # for isomers.xyz
    x = {"file": "tests/data/isomers.xyz",
         "bond_broken": [3, 2], "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    Result = Sts_topo.method_broken_bond_H(args)
    assert Result == [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]


def test_topo_Broken_bond():
    # for crest_conformers.xyz
    x = {"file": "tests/data/crest_conformers.xyz",
         "bond_broken": [52, 55], "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    Result = Sts_topo.method_broken_bond(args)
    assert Result == [1, 2, 3, 6, 7, 8, 13, 15, 16, 19, 21,
                      23, 24, 27, 30, 33, 35, 37, 40, 44, 45, 50, 51, 52]
    # for isomers.xyz
    x = {"file": "tests/data/isomers.xyz",
         "bond_broken": [3, 2], "print": True, "debug": False}
    args = argparse.Namespace(**x)
    Sts_topo: Topo = Topo(Path(args.file))
    Result = Sts_topo.method_broken_bond(args)
    assert Result == [1, 3, 6, 9, 10, 17]
