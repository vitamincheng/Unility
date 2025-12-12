from pathlib import Path
from ase import Atoms
import numpy.typing as npt
import numpy as np
import pytest
from icecream import ic
from censo_ext.Tools.ml4nmr import read_mol_neighbors, read_mol_neighbors_bond_order
from censo_ext.Tools.utility import IsExists_DirFileName
from censo_ext.Tools.xyzfile import GeometryXYZs

mol: Atoms | list[Atoms]
idx_neighbors: dict[int, npt.NDArray[np.int64]]


def test_m4nmr_miss_args():
    FileName: Path = Path("tests/data/crest_conformers0000.xyz")
    with pytest.raises(SystemExit) as e:
        IsExists_DirFileName(FileName)
        xyzFile: GeometryXYZs = GeometryXYZs(FileName)
        xyzFile.method_read_xyz()
        mol, idx_neighbors = read_mol_neighbors(
            xyzFile)
    assert e.type is SystemExit
    assert e.value.code == 0

    with pytest.raises(SystemExit) as e:
        xyzFile: GeometryXYZs = GeometryXYZs(FileName)
        xyzFile.method_read_xyz()
        mol, idx_neighbors, bond_order = read_mol_neighbors_bond_order(
            xyzFile)
    assert e.type is SystemExit
    assert e.value.code == 0


def test_m4nmr_read_mol_neighbors():
    # for crest_conformers.xyz
    inFile = Path("tests/data/crest_conformers.xyz")

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    mol, idx_neighbors = read_mol_neighbors(
        xyzFile)
    xyzFile: GeometryXYZs = GeometryXYZs(
        Path("tests/data/crest_conformers.xyz"))
    xyzFile.method_read_xyz()
    mol, idx_neighbors, bond_order = read_mol_neighbors_bond_order(
        xyzFile)
    nAtoms: int = (len([atom for atom in mol]))
    assert nAtoms == 73
    assert [atom.number for atom in mol] == [6, 6, 6, 1, 1, 6, 6, 6, 1, 1, 1, 1, 8, 1, 6, 6, 1, 1, 6, 1, 6, 1, 6, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1,
                                             6, 1, 6, 1, 1, 6, 1, 1, 1, 6, 6, 1, 1, 1, 1, 6, 6, 6, 1, 1, 6, 6, 6, 1, 1, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    _, len_neighbors, bond_order = (mol, len(idx_neighbors), bond_order)
    assert len_neighbors == 73
    assert bond_order == {1: 0, 2: 2, 3: 1, 6: 0, 7: 2, 8: 2, 15: 1, 16: 2, 19: 1, 21: 0, 23: 2, 24: 2, 27: 2,
                          30: 0, 33: 1, 35: 2, 37: 2, 40: 1, 44: 1, 45: 3, 50: 3, 51: 1, 52: 1, 55: 1, 56: 3, 57: 1, 60: 3, 61: 3}

    # for isomers.xyz
    inFile = Path("tests/data/isomers.xyz")

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    mol, idx_neighbors = read_mol_neighbors(xyzFile)
    ic(mol, idx_neighbors)
    xyzFile: GeometryXYZs = GeometryXYZs(Path("tests/data/isomers.xyz"))
    xyzFile.method_read_xyz()
    mol, idx_neighbors, bond_order = read_mol_neighbors_bond_order(
        xyzFile)
    nAtoms = (len([atom for atom in mol]))
    assert nAtoms == 17
    assert [atom.number for atom in mol] == [
        6, 6, 6, 1, 1, 6, 1, 1, 6, 6, 1, 1, 1, 1, 1, 1, 8]

    _, len_neighbors, bond_order = (mol, len(idx_neighbors), bond_order)
    assert len_neighbors == 17
    assert bond_order == {1: 2, 2: 2, 3: 2, 6: 0, 9: 2, 10: 2}
