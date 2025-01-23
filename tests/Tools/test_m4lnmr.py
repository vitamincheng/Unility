import pathlib
import argparse
from pathlib import Path
from tkinter import W
import pytest
from icecream import ic
from Tools.ml4nmr import *


def test_m4nmr():

    mol, neighbors = read_mol_neighbors(
        Path("tests/data/crest_conformers.xyz"))
    ic(mol, neighbors)
    mol, neighbors, bond_order = read_mol_neighbors_bond_order(
        Path("tests/data/crest_conformers.xyz"))
    atoms_nums = (len([atom for atom in mol]))
    assert atoms_nums == 73
    atoms_list = ([atom.number for atom in mol])  # type: ignore
    assert atoms_list == [6, 6, 6, 1, 1, 6, 6, 6, 1, 1, 1, 1, 8, 1, 6, 6, 1, 1, 6, 1, 6, 1, 6, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1,
                          6, 1, 6, 1, 1, 6, 1, 1, 1, 6, 6, 1, 1, 1, 1, 6, 6, 6, 1, 1, 6, 6, 6, 1, 1, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    _, len_neighbors, bond_order = (mol, len(neighbors), bond_order)
    assert len_neighbors == 73
    assert bond_order == {1: 0, 2: 2, 3: 1, 6: 0, 7: 2, 8: 2, 15: 1, 16: 2, 19: 1, 21: 0, 23: 2, 24: 2, 27: 2,
                          30: 0, 33: 1, 35: 2, 37: 2, 40: 1, 44: 1, 45: 3, 50: 3, 51: 1, 52: 1, 55: 1, 56: 3, 57: 1, 60: 3, 61: 3}
