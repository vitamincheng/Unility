#!/usr/bin/env python3
# |                                       [08.18.2024] vitamin.cheng@gmail.com
#  Need ASE Library
########## GLOBAL DECLARATIONS ##########

# use slightly modified covalent radii from ase for neighbor recognition
from pathlib import Path
from ase.atoms import Atoms
import numpy as np
import numpy.typing as npt
from icecream import ic
from censo_ext.Tools.utility import IsExist
import os
from ase.data import covalent_radii
custom_radii: npt.NDArray = covalent_radii.copy()
custom_radii[3] -= 0.15   # reduce radius of Li
custom_radii[6] -= 0.05   # reduce radius of C

# Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197).
# Values for metals decreased by 10 %.
# This was copied from gituhb project dftd3/tad-dftd3/src/tad_dftd3/data.py
covalent_rad_2009: npt.NDArray = np.array([
    0.00,                                                # None
    0.32, 0.46,                                           # H,He
    1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,             # Li-Ne
    1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,             # Na-Ar
    1.76, 1.54,                                           # K,Ca
    1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,   # Sc-Zn
    1.12, 1.09, 1.15, 1.10, 1.14, 1.17,                       # Ga-Kr
    1.89, 1.67,                                           # Rb,Sr
    1.47, 1.39, 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23,   # Y-Cd
    1.28, 1.26, 1.26, 1.23, 1.32, 1.31,                       # In-Xe
    2.09, 1.76, 1.62,                                      # Cs-La
    1.47, 1.58, 1.57, 1.56, 1.55, 1.51, 1.52,                  # Ce-Gd
    1.51, 1.50, 1.49, 1.49, 1.48, 1.53, 1.46,                  # Tb-Lu
    1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,        # Hf-Hg
    1.30, 1.30, 1.36, 1.31, 1.38, 1.42,                       # Tl-Rn
    2.01, 1.81, 1.67,                                      # Fr-Ac
    1.58, 1.52, 1.53, 1.54, 1.55, 1.49, 1.49,                  # Th-Cm
    1.51, 1.51, 1.48, 1.50, 1.56, 1.58, 1.45,                  # Bk-Lr
    1.41, 1.34, 1.29, 1.27, 1.21, 1.16, 1.15, 1.09, 1.22,        # Rf-Cn
    1.36, 1.43, 1.46, 1.58, 1.48, 1.57                        # Nh-Og
])

# D3 covalent radii used to construct the coordianation number
covalent_rad_d3 = 4.0 / 3.0 * covalent_rad_2009

########## END GLOBAL DECLARATIONS ##########


# def get_number_nmr(outlist: list) -> int:
#    """Get number of calculated NMR nuclei (not necessarily only 1H and 13C).
#
#    outlist: list of lines read from a calculation output file
#    """
#
#    n_nmr: int = 0
#    for line in outlist[outlist.index('                            *     ORCA property calculations      *\n'):]:
#        if 'Number of nuclei for epr/nmr' in line:
#            n_nmr = int(line.split()[-1])
#    if n_nmr == 0:
#        print("ERROR: Number of calculated NMR nuclei was not found in ORCA output!")
#        exit()
#
#    return n_nmr
#
#
# def read_orca(filename: Path = Path("orcaS.out")) -> tuple:
#    """Read the DFT output from ORCA calculations."""
#
#    data: list[str] = open(filename, 'r').readlines()
#
#    # get number of calculated NMR nuclei
#    n_nmr: int = get_number_nmr(data)
#
#    # get the range in which the shieldings are listed
#    start: int = data.index('CHEMICAL SHIELDING SUMMARY (ppm)\n') + 6
#    end: int = start + n_nmr
#
#    shieldings: list = []
#    for line in data[start:end]:
#        tmp: list[str] = line.split()
#        shieldings.append({
#            # ORCA starts at 0 counting the nuc numbers
#            'idx_nuclei': int(tmp[0]) + 1,
#            'element': tmp[1].upper(),
#            'Shielding': float(tmp[2])
#        })
#
#    # make sure the shieldings are ordered according to the atom numbering (store as tuple)
#    return tuple(sorted(shieldings, key=lambda s: s['idx_nuclei']))
#
#
# def getref_orca(filename: Path = Path("orcaS.out")) -> tuple[list, list]:
#    """Get the reference shieldings from the DFT output (ORCA).
#
#    ATTENTION: all Hs and all Cs are averaged, works e.g. for TMS and CH4.
#    """
#    data: list[str] = open(filename, 'r').readlines()
#
#    # get number of calculated NMR nuclei
#    n_nmr: int = get_number_nmr(data)
#
#    # get the range in which the shieldings are listed
#    start: int = data.index('CHEMICAL SHIELDING SUMMARY (ppm)\n') + 6
#    end: int = start + n_nmr
#
#    list_h: list = []
#    list_c: list = []
#    idx_h, idx_c = 0, 0
#
#    for line in data[start:end]:
#        tmp: list = line.split()
#        if tmp[1].upper() == 'H':
#            list_h.append(float(tmp[2]))
#            idx_h += 1
#        if tmp[1].upper() == 'C':
#            list_c.append(float(tmp[2]))
#            idx_c += 1
#
#    return list_h, list_c


def read_mol_neighbors(DirFileName: Path):
    """Read the molecule and return mol (ase.Atoms object) and dict neighbors."""

    # read the .xyz coordinates from the molecular structures
    import ase.io
    from ase import neighborlist
    IsExist(DirFileName)
    mol: Atoms | list[Atoms] = ase.io.read(str(DirFileName), format='xyz')

    # use covalent radii as thresholds for neighbor determination (what about vdW radii?)
    cutoffs = [custom_radii[atom.number] for atom in mol]  # type: ignore

    # build neighbor list and write list of neighboring atoms to the dict neighbors
    nl = neighborlist.build_neighbor_list(
        mol, cutoffs, self_interaction=False, bothways=True)
    neighbors: dict[int, np.ndarray] = {}

    for idx in range(len(mol)):
        # nl.get_neighbors(i) returns [0]: indices and [1]: offsets
        indices = nl.get_neighbors(idx)[0]
        # add 1 to key and to value to start counting of atoms at 1
        neighbors[idx+1] = indices+1

        # exit if an H atom has not exactly 1 neighbor
        if mol.get_atomic_numbers()[idx] == 1 and len(neighbors[idx+1]) != 1:  # type: ignore # nopep8
            print(f"ERROR: H atom {idx+1} has not exactly one neighbor! File in: {DirFileName}")  # nopep8
            exit(0)

    return mol, neighbors


def read_mol_neighbors_bond_order(DirfileName: Path = Path("crest_conformers.xyz")):
    """Read the molecule and return mol (ase.Atoms object) and dict neighbors."""

    # read the .xyz coordinates from the molecular structures
    mol, neighbors = read_mol_neighbors(DirfileName)

    idx_H_atom: list[int] = [idx+1 for idx, i in enumerate(mol) if i.symbol == "H"]  # type: ignore # nopep8
    idx_C_atom: list[int] = [idx+1 for idx, i in enumerate(mol) if i.symbol == "C"]  # type: ignore # nopep8
    bond_order: dict[int, int] = {}
    for idx in neighbors.keys():
        count = 0
        for idy in neighbors[idx]:
            if idy in idx_H_atom:
                count = count + 1
        if idx in idx_C_atom:
            bond_order[idx] = count

    return mol, neighbors, bond_order
