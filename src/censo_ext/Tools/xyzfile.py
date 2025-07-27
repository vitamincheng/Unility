#!/usr/bin/env python3

#  Modification             [08.26.2024] vitamin.cheng@gmail.com
#

from __future__ import annotations
from pathlib import Path
from typing import Self
import os
import sys
import numpy as np
import numpy.typing as npt
import argparse
from icecream import ic
import copy


class Geometry():
    '''Stores all of the data in an xyz file (single confromer) '''

    def __init__(self, names: dict[int, str], coord: list[npt.NDArray], extras: list[list[str]], comment: str = "", energy: float = 0, ncluster: int = 0) -> None:

        self.names: dict[int, str] = names          # atom's name   H Li Na K B C O S F Cl # nopep8
        self.coord: list[npt.NDArray] = coord       # coordinates of every atom
        self.nAtoms: int = len(names)               # numbers of atom
        self.comment: str = comment                 # Energy =   Eh   #Cluster  :i         # nopep8
        self.comment_energy: float = energy         # Energy (Eh)
        self.comment_ncluster: int = ncluster       # index of Clusters
        self.mass: npt.NDArray
        self.extras: list[list[str]] = extras
        self.com: npt.NDArray
        self.inertia: npt.NDArray

    def __add__(self, other: Self) -> Geometry:
        import copy
        result: Geometry = copy.deepcopy(self)
        for key in other.names.keys():
            result.names[len(self.names) + key] = other.names[key]
        # coordinates of every atom
        result.coord = self.coord + other.coord
        result.extras = self.extras + other.extras
        result.nAtoms = len(self.names) + len(other.names)
        return result

    def __repr__(self) -> str:
        result: str = ""
        result += f"{self.nAtoms}\n{self.comment}\n"
        # ic(self.extras, len(self.extras))
        for i in range(self.nAtoms):
            extra: str = "   ".join(
                self.extras[i]) if self.extras[i] != [] else ""
            result += f'{self.names[i+1]:>3s}    {self.coord[i][0]: 14.10f}'  # nopep8
            result += f'    {self.coord[i][1]: 14.10f}    {self.coord[i][2]: 14.10f}    {extra}\n'  # nopep8
        return result

    def get_comment_energy(self) -> float:
        return self.comment_energy

    def method_translate_xyz(self, delta: npt.NDArray) -> Geometry:
        result: Geometry = copy.deepcopy(self)
        for x in result.coord:
            x += delta
        return result

    def method_molecules_separation_xyz(self, list_idx: list) -> bool:
        new_names: dict = {}
        idx_names: int = 1
        for name in self.names.copy():
            if name in list_idx:
                new_names[idx_names] = self.names[name]
                idx_names = idx_names + 1
            else:
                del self.names[name]
        self.names = new_names
        self.nAtoms = len(self.names)
        idx0_St: list = [x-1 for x in list_idx]
        self.coord = (np.array(self.coord)[idx0_St]).tolist()
        return True

    def method_idx_molecules_xyz(self, fileName: Path) -> list[int]:
        from censo_ext.Tools.topo import Topo
        import censo_ext.Tools.ml4nmr as ml4nmr
        mol, neighbors = ml4nmr.read_mol_neighbors(fileName)
        # ic(neighbors)
        idx_molecule: set[int] = set([*range(1, len(neighbors)+1, 1)])
        molecules: list = []
        while (len(idx_molecule) != 0):
            idx_atoms = 0
            Hydrogen: list[int] = []
            for x in idx_molecule:
                if len(neighbors[x]) == 1 or len(neighbors[x]) == 0:
                    idx_atoms = x
                    break
            if len(neighbors[idx_atoms]) == 1:
                x = {"file": fileName, "bonding": int(
                    idx_atoms), "print": False, "debug": False}
                Sts_topo: Topo = Topo(x["file"])  # type: ignore
                neighbors_bonding: list[int] = Sts_topo.method_bonding(
                    argparse.Namespace(**x))
                x = {"file": fileName, "bond_broken": [
                    neighbors_bonding[0], idx_atoms], "print": False, "debug": False}
                Hydrogen = Sts_topo.method_broken_bond_H(
                    argparse.Namespace(**x))
            elif len(neighbors[idx_atoms]) == 0:
                Hydrogen = [idx_atoms]
            molecules.append(Hydrogen)
            idx_molecule = idx_molecule.difference(set(Hydrogen))
        return molecules

    def method_computeCOM(self) -> npt.NDArray:
        '''
        Returns the center of mass of the geometry.
        '''
        self.method_update_masses()
        if (len(self.com) == 3):
            return self.com
        else:
            self.com = np.dot(self.mass, self.coord) / np.sum(self.mass)
            return self.com

    def method_computeInertia(self) -> npt.NDArray:
        '''
        Returns the moment of inertia tensor
        '''
        com = self.method_computeCOM()
        data = self.coord - com

        self.inertia = -np.einsum("ax,a,ay->xy", data, self.mass, data)
        return self.inertia

    def method_update_masses(self) -> None:
        self.mass: npt.NDArray = np.array([])
        for idx, x in self.names.items():
            from censo_ext.Tools.Parameter import masses
            self.mass = np.append(self.mass, masses[x.lower()], axis=None)

    def method_rewrite_comment(self) -> None:
        self.comment = " Energy = "+" "*7 + \
            f"{self.comment_energy:.10f}"+" Eh"+" "*8 + \
            "#Cluster:     "+str(self.comment_ncluster)

    def method_update_comment(self) -> None:
        '''
        Analysis of comment
        '''
        comments: list = self.comment.replace("a.u.", "").replace("Eh", "").replace("Energy=", "").replace("Energy =", "").replace(
            "Energy  =", "").replace("energy:", "").replace("Energy:", "").split()
        if comments == []:
            print(" Your xyz file have not any about Energy and Cluster !!!")
            print(" We will set Energy = 0 in your xyz file")
            self.comment_energy, self.comment_ncluster = 0, 0
            self.method_rewrite_comment()
            return

        from censo_ext.Tools.utility import function_is_float
        if function_is_float(comments[0]):
            if len(comments) >= 3:
                if comments[1] == "#Cluster:" and function_is_float(comments[2]):
                    self.comment_energy, self.comment_ncluster = float(comments[0]), int(comments[2])  # nopep8
                else:
                    self.comment_energy, self.comment_ncluster = float(comments[0]), 0  # nopep8

                # if comments[4] == "E" and function_is_float(comments[5]):
                #    self.comment_energy, self.comment_ncluster = float(comments[5]), 0  # nopep8

            else:
                self.comment_energy, self.comment_ncluster = float(comments[0]), 0  # nopep8
        else:
            print("something wrong in your xyz file !!!")
            print(comments)
            print(" Quit the program !!!")
            ic()
            exit(1)

        self.method_rewrite_comment()
        return

    def method_comment_new(self, idx1: int) -> None:
        self.comment_ncluster = idx1
        self.method_rewrite_comment()


class GeometryXYZs():
    '''
    stores all the data in an xyz file (multi-conformers)
    '''

    def __init__(self, fileName: Path = Path("")) -> None:
        self.__filename: Path = Path(fileName)
        self.Sts: list[Geometry] = list()

    def __len__(self) -> int:
        return int(len(self.Sts))

    def set_filename(self, fileName: Path) -> None:
        self.__filename = Path(fileName)

    def method_translate_cut_xyzs(self, delta: npt.NDArray, cut: int) -> GeometryXYZs:
        result: GeometryXYZs = GeometryXYZs()
        if len(self) == 1:
            for x in np.linspace(0, 1, num=cut, endpoint=True):
                dSt = self.Sts[0].method_translate_xyz(delta=delta*x)
                result.Sts.append(copy.deepcopy(dSt))
            return result
        else:
            print(" More than 1 in your xyz file ")
            print(" All xyzs will to save your xyz file")
            for St in (self.Sts):
                for x in np.linspace(0, 1, num=cut, endpoint=True):
                    dSt: Geometry = St.method_translate_xyz(delta=delta*x)
                    result.Sts.append(copy.deepcopy(dSt))
            return result

    def method_translate_xyzs(self, delta: npt.NDArray) -> GeometryXYZs:
        result: GeometryXYZs = GeometryXYZs()
        for x in self.Sts:
            x: Geometry = x.method_translate_xyz(delta=delta)
            result.Sts.append(copy.deepcopy(x))
        return result

    def __add__(self, var2: Self) -> GeometryXYZs:

        var3: GeometryXYZs
        if len(var2) == 1:
            var3: GeometryXYZs = copy.deepcopy(self)
            var3.Sts.append(var2.Sts[0])

            var4: GeometryXYZs = GeometryXYZs()
            for idx in range(0, len(self), 1):
                var4.Sts.append(
                    var3.Sts[idx] + var3.Sts[-1])
            return var4
        else:
            print("Too much xyzs structures in your xyz file")
            ic()
            exit(1)

    def method_idx_molecules_xyzs(self, idx1: int = 1) -> bool:
        fileName: Path = Path("~temp.xyz")
        self.set_filename(fileName)
        self.method_save_xyz([idx1])
        list_idx: list = self.Sts[idx1 -
                                  1].method_idx_molecules_xyz(fileName)
        from censo_ext.Tools.utility import delete_all_files
        delete_all_files(fileName)
        for x in list_idx:
            self.Sts.append(copy.deepcopy(self.Sts[idx1-1]))
        # ic(len(self.structures))
        for idx, x in enumerate(list_idx):
            if self.Sts[idx1+idx].method_molecules_separation_xyz(x):
                pass
            else:
                print(" Something wrong in your Molecule Separation xyz file")
                print(" Exit to the program !!!")
                ic()
                exit(1)

        path = "Separation"
        isExist: bool = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        for idx0, x in enumerate(list_idx):
            self.set_filename(Path(path) / Path(str(idx0+1)+".xyz"))
            self.method_save_xyz([idx1+idx0+1])
        return True

    def method_Sts_extend(self, Sts_: list[Geometry]) -> None:
        raise NotImplementedError
        # self.Sts.extend(Sts_)

    def method_Sts_append(self, St: Geometry) -> None:
        raise NotImplementedError
        # self.Sts.append(St)

    def method_read_xyz(self) -> None:
        '''reads xyz file and save the data into self geometry object'''
        from censo_ext.Tools.utility import IsExist
        IsExist(self.__filename)

        with open(self.__filename, "r") as f:
            line: str = f.readline()
            while line != "":
                nAtoms = int(line)
                comment: str = f.readline().rstrip()
                names: dict[int, str] = dict()
                coords: list[npt.NDArray] = list()
                extras: list[list[str]] = list()

                for i in range(nAtoms):
                    line = f.readline()
                    data: list[str] = line.split()
                    name, x, y, z = data[0:4]
                    extra: list[str] = data[4:]

                    names[i+1] = name.capitalize()
                    coords.append(np.array([float(x), float(y), float(z)]))
                    if extra:
                        extras.append(extra)
                    else:
                        extras.append(list([]))

                self.Sts.append(Geometry(
                    names, coords, comment=comment, extras=extras))
                line = f.readline()
        self.method_comment_keep()

    def method_save_xyz(self, idx1_list: list[int]) -> None:
        original_stdout = sys.stdout
        with open(self.__filename, "w") as f:
            sys.stdout = f
            self.method_print(idx1_list)
        sys.stdout = original_stdout

    def method_save_xyz_append(self, idx1_list: list) -> None:  # append to old xyz file
        original_stdout = sys.stdout
        with open(self.__filename, "a") as f:
            sys.stdout = f
            self.method_print(idx1_list)
        sys.stdout = original_stdout

    def method_print(self, idx1_St: list[int]) -> None:
        if (idx1_St == []):
            idx0_St: list = [*range(len(self))]
        else:
            idx0_St: list = [x-1 for x in idx1_St]
        for key in idx0_St:
            print(self.Sts[key], end="")

    def method_comment_keep(self) -> None:
        for St in self.Sts:
            if St.method_update_comment():
                St.method_update_comment()

    def method_comment_new(self) -> None:
        for idx0, St in enumerate(self.Sts):
            if St.method_update_comment():
                St.method_update_comment()
            St.method_comment_new(idx0+1)

    def method_rewrite_comment(self) -> None:
        for St in self.Sts:
            if St.method_rewrite_comment():
                St.method_rewrite_comment()

    def get_comment_energy(self) -> list:
        result: list = []
        for St in self.Sts:
            if St.get_comment_energy():
                result.append(St.get_comment_energy())
        return result

    def method_ensoGenFlexible(self, args, thermo_list):
        import numpy.lib.recfunctions as rfn
        from censo_ext.Tools.Parameter import Eh, FACTOR
        # dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'),
        #       ('Energy', '<f8'), ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8')]

        # Column 8 is Total Gibbs Free Energy (Eh) = Energy + mRRHO
        # ic(thermo_list)
        TEMP: float = args.temp
        np_enso: npt.NDArray = np.zeros((len(self.Sts),), dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'),
                                                                 ('Energy', '<f8'), ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8')])
        np_enso['ONOFF'] = 1
        np_enso['gi'] = 1.000
        np_enso['Gsolv'] = 0.00000000
        np_enso['NMR'] = np.arange(1, len(self.Sts)+1)
        np_enso['CONF'] = np.arange(1, len(self.Sts)+1)
        np_enso['mRRHO'] = np.array(thermo_list)
        np_enso['Energy'] = np.array(
            [a.comment_energy for a in self.Sts], dtype=[('Energy', 'f8')])
        Total: npt.NDArray = np.array(
            np_enso['Energy']+np_enso['mRRHO'], dtype=[('Total', 'f8')])
        np_enso = rfn.merge_arrays((np_enso, Total), flatten=True)

        # Gibbs_min is lowest energy of Gibbs Free Energy
        Gibbs_min: np.float64 = np_enso['Total'].min()

        # Column 9 is delta Gibbs Free Energy (kcal/mol)
        Gibbs: npt.NDArray = np.array(
            (np_enso['Total']-Gibbs_min)*Eh, dtype=[('Gibbs', 'f8')])
        np_enso = rfn.merge_arrays((np_enso, Gibbs), flatten=True)

        # Column 1o is Qi (each CONFS)
        Qi: npt.NDArray = np.array(
            np.exp(-np_enso['Gibbs']/(TEMP*FACTOR)), dtype=[('Qi', 'f8')])
        np_enso = rfn.merge_arrays((np_enso, Qi), flatten=True)

        # Qall is sum of Qi
        Qall: np.float64 = np.sum(np_enso['Qi'])

        # Column 11 is percentage of each CONFS
        NEW_BW: npt.NDArray = np.array(
            np_enso['Qi']/Qall, dtype=[('NEW_BW', 'f8')])
        np_enso = rfn.merge_arrays((np_enso, NEW_BW), flatten=True)
        # copy BW and delete residual parameter
        np_enso['BW'] = np_enso['NEW_BW']
        names_anmr: list = list()
        if np_enso.dtype.names != None:
            names_anmr = list(np_enso.dtype.names)
        return np_enso[names_anmr[:8]]
