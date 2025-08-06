#!/usr/bin/env python

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

    def __init__(self, names: dict[int, str], coord: list[npt.NDArray[np.float64]], extras: list[list[str]], comment: str = "", energy: float = 0, nClusters: int = 0) -> None:
        """ 
        Initialize a Geometry object with atom names, coordinates, and metadata.

        Args:
            names (dict[int, str]): Dictionary mapping atom indices to names.
            coord (list[npt.NDArray[np.float64]]): List of atomic coordinates.
            extras (list[list[str]]): Extra data for each atom (e.g., charges).
            comment (str): Optional comment string.
            energy (float): Optional energy value (in Eh).
            nClusters (int): Optional cluster index.
        """

        self.names: dict[int, str] = names                  # atom's name   H Li Na K B C O S F Cl # nopep8
        self.coord: list[npt.NDArray[np.float64]] = coord
        # coordinates of every atom #type:ignore #nopep8
        self.nAtoms: int = len(names)                       # numbers of atom
        self.comment: str = comment                         # Energy =   Eh   #Cluster  :i         # nopep8
        self.comment_energy: float = energy                 # Energy (Eh)
        self.comment_nClusters: int = nClusters             # index of Clusters
        self.mass: npt.NDArray
        self.extras: list[list[str]] = extras
        self.com: npt.NDArray
        self.inertia: npt.NDArray

    def __add__(self, other: Self) -> Geometry:
        """ 
        Concatenate two Geometry objects by combining atoms and coordinates.

        Returns:
            Geometry: New Geometry instance with combined data.
        """

        import copy
        res: Geometry = copy.deepcopy(self)
        for key in other.names.keys():
            res.names[len(self.names) + key] = other.names[key]
        # coordinates of every atom
        res.coord = self.coord + other.coord
        res.extras = self.extras + other.extras
        res.nAtoms = len(self.names) + len(other.names)
        return res

    def __repr__(self) -> str:
        """ 
        Generate a string representation of the Geometry object in XYZ format.

        Returns:
            str: Formatted XYZ string with atom names, coordinates, and extras.
        """

        res: str = ""
        res += f"{self.nAtoms}\n{self.comment}\n"
        for idx0 in range(self.nAtoms):
            extra: str = "   ".join(
                self.extras[idx0]) if self.extras[idx0] != [] else ""
            res += f'{self.names[idx0+1]:>3s}    {self.coord[idx0][0]: 14.10f}'  # nopep8
            res += f'    {self.coord[idx0][1]: 14.10f}    {self.coord[idx0][2]: 14.10f}    {extra}\n'  # nopep8
        return res

    def get_comment_energy(self) -> float:
        """
        Retrieve the energy value from the comment field.

        Returns:
            float: Energy value in Eh units.
        """

        return self.comment_energy

    def method_translate_xyz(self, delta: npt.NDArray) -> Geometry:
        """
        Translate all atoms by a given vector.

        Args:
            delta (npt.NDArray): Translation vector (x, y, z).

        Returns:
            Geometry: Translated Geometry instance.
        """

        res: Geometry = copy.deepcopy(self)
        for x in res.coord:
            x += delta
        return res

    def method_molecules_separation_xyz(self, idx1_Select_Names: list[int]) -> bool:
        """
        Filter atoms by selected indices and update the Geometry.

        Args:
            idx1_Select_Names (list[int]): List of atom indices to retain.

        Returns:
            bool: True if operation succeeded.
        """

        new_names: dict[int, str] = {}
        x: int = 1
        for names_key in self.names.copy():
            if names_key in idx1_Select_Names:
                new_names[x] = self.names[names_key]
                x = x + 1
            else:
                del self.names[names_key]
        self.names = new_names
        self.nAtoms = len(self.names)
        idx0_St: list[int] = [x-1 for x in idx1_Select_Names]
        self.coord = (np.array(self.coord)[idx0_St]).tolist()
        return True

    def method_idx_molecules_xyz(self, fileName: Path) -> list[list[int]]:
        """
        Identify molecular clusters from a topology file.

        Args:
            fileName (Path): Path to topology file.

        Returns:
            list[list[int]]: List of molecule atom index lists.
        """

        from censo_ext.Tools.topo import Topo
        import censo_ext.Tools.ml4nmr as ml4nmr
        neighbors: dict[int, npt.NDArray[np.int64]]
        _, neighbors = ml4nmr.read_mol_neighbors(fileName)
        idx1_molecule: set[int] = set([*range(1, len(neighbors)+1, 1)])
        molecules: list[list[int]] = []
        while (len(idx1_molecule) != 0):
            idx_atoms: int = 0
            Hydrogen: list[int] = []
            for x in idx1_molecule:
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
            idx1_molecule = idx1_molecule.difference(set(Hydrogen))
        return molecules

    def method_computeCOM(self) -> npt.NDArray:
        """
        Calculate the center of mass (COM) of the molecule.

        Returns:
            npt.NDArray: COM coordinates.
        """

        self.method_update_masses()
        if (len(self.com) == 3):
            return self.com
        else:
            self.com = np.dot(self.mass, self.coord) / np.sum(self.mass)
            return self.com

    def method_computeInertia(self) -> npt.NDArray:
        """
        Calculate the moment of inertia tensor.

        Returns:
            npt.NDArray: 3x3 inertia tensor.
        """

        data = self.coord - self.method_computeCOM()

        self.inertia = -np.einsum("ax,a,ay->xy", data, self.mass, data)
        return self.inertia

    def method_update_masses(self) -> None:
        """
        Update atomic masses based on element names.
        """

        self.mass = np.array([])
        for _, x in self.names.items():
            from censo_ext.Tools.Parameter import masses
            self.mass = np.append(self.mass, masses[x.lower()], axis=None)

    def method_rewrite_comment(self) -> None:
        """
        Format and update the comment field with energy and cluster info.
        """

        self.comment = " Energy = "+" "*7 + \
            f"{self.comment_energy:.10f}"+" Eh"+" "*8 + \
            "#Cluster:     "+str(self.comment_nClusters)

    def method_update_comment(self) -> None:
        """
        Parse and extract energy and cluster information from the comment field.
        """

        comments: list[str] = self.comment.replace("a.u.", "").replace("Eh", "").replace("Energy=", "").replace("Energy =", "").replace(
            "Energy  =", "").replace("energy:", "").replace("Energy:", "").split()
        if comments == []:
            print(" Your xyz file have not any about Energy and Cluster !!!")
            print(" We will set Energy = 0 in your xyz file")
            self.comment_energy, self.comment_nClusters = 0, 0
            self.method_rewrite_comment()
            return

        from censo_ext.Tools.utility import function_is_float
        if function_is_float(comments[0]):
            if len(comments) >= 3:
                if comments[1] == "#Cluster:" and function_is_float(comments[2]):
                    self.comment_energy, self.comment_nClusters = float(comments[0]), int(comments[2])  # nopep8
                else:
                    self.comment_energy, self.comment_nClusters = float(comments[0]), 0  # nopep8

                # if comments[4] == "E" and function_is_float(comments[5]):
                #    self.comment_energy, self.comment_nClusters = float(comments[5]), 0  # nopep8

            else:
                self.comment_energy, self.comment_nClusters = float(comments[0]), 0  # nopep8
        else:
            print("something wrong in your xyz file !!!")
            print(comments)
            print(" Quit the program !!!")
            ic()
            raise ValueError(" something wrong in your xyz file !!! ")

        self.method_rewrite_comment()
        return

    def method_comment_new(self, idx1: int) -> None:
        """
        Set a new cluster index in the comment field.

        Args:
            idx1 (int): Cluster index to set.
        """

        self.comment_nClusters = idx1
        self.method_rewrite_comment()


class GeometryXYZs():

    def __init__(self, fileName: Path = Path("")) -> None:
        """
        Initialize a GeometryXYZs object to manage multiple Geometry instances.

        Args:
            fileName (Path): Path to XYZ file.
        """

        self.__filename: Path = Path(fileName)
        self.Sts: list[Geometry] = list()

    def __len__(self) -> int:
        return int(len(self.Sts))

    def set_filename(self, fileName: Path) -> None:
        self.__filename = Path(fileName)

    def method_translate_cut_xyzs(self, delta: npt.NDArray[np.float64], cut: int) -> GeometryXYZs:
        """
        Generate interpolated GeometryXYZs by translating along a vector.

        Args:
            delta (npt.NDArray[np.float64]): Translation vector.
            cut (int): Number of interpolation points.

        Returns:
            GeometryXYZs: Interpolated structures.
        """

        res: GeometryXYZs = GeometryXYZs()
        if len(self) == 1:
            for x in np.linspace(0, 1, num=cut, endpoint=True):
                dSt = self.Sts[0].method_translate_xyz(delta=delta*x)
                res.Sts.append(copy.deepcopy(dSt))
            return res
        else:
            print(" More than 1 in your xyz file ")
            print(" All xyzs will to save your xyz file")
            for St in self.Sts:
                for x in np.linspace(0, 1, num=cut, endpoint=True):
                    dSt: Geometry = St.method_translate_xyz(delta=delta*x)
                    res.Sts.append(copy.deepcopy(dSt))
            return res

    def method_translate_xyzs(self, delta: npt.NDArray[np.float64]) -> GeometryXYZs:
        """
        Translate all Geometry instances in the collection.

        Args:
            delta (npt.NDArray[np.float64]): Translation vector.

        Returns:
            GeometryXYZs: Translated structures.
        """

        Res: GeometryXYZs = GeometryXYZs()
        for x in self.Sts:
            x: Geometry = x.method_translate_xyz(delta=delta)
            Res.Sts.append(copy.deepcopy(x))
        return Res

    def __add__(self, Var: Self) -> GeometryXYZs:
        """
        Concatenate two GeometryXYZs objects.

        Args:
            Var (GeometryXYZs): Other GeometryXYZs instance.

        Returns:
            GeometryXYZs: Concatenated structures.
        """

        GeoXYZs: GeometryXYZs
        if len(Var) == 1:
            GeoXYZs = copy.deepcopy(self)
            GeoXYZs.Sts.append(Var.Sts[0])

            res: GeometryXYZs = GeometryXYZs()
            for idx in range(0, len(self), 1):
                res.Sts.append(
                    GeoXYZs.Sts[idx] + GeoXYZs.Sts[-1])
            return res
        else:
            print("Too much xyzs structures in your xyz file")
            ic()
            raise ValueError("Too much xyzs structures in your xyz file")

    def method_idx_molecules_xyzs(self, idx1: int = 1) -> bool:
        """
        Split molecules in all Geometry instances using topology data.

        Args:
            idx1 (int): Starting index for molecule separation.

        Returns:
            bool: True if operation succeeded.
        """

        fileName: Path = Path("~temp.xyz")
        self.set_filename(fileName)
        self.method_save_xyz([idx1])
        list_idx: list[list[int]] = self.Sts[idx1 -
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
                raise ValueError(
                    " Something wrong in your Molecule Separation xyz file")

        path = "Separation"
        isExist: bool = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        for idx0, x in enumerate(list_idx):
            self.set_filename(Path(path) / Path(str(idx0+1)+".xyz"))
            self.method_save_xyz([idx1+idx0+1])
        return True

    def method_Sts_extend(self, Sts_: list[Geometry]) -> None:
        raise NotImplementedError("Under Construct")
        # self.Sts.extend(Sts_)

    def method_Sts_append(self, St: Geometry) -> None:
        raise NotImplementedError("Under Construct")
        # self.Sts.append(St)

    def method_read_xyz(self) -> None:
        """
        Read XYZ file and populate the GeometryXYZs collection.
        """

        from censo_ext.Tools.utility import IsExists_DirFileName
        IsExists_DirFileName(self.__filename)

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
        """
        Save selected Geometry instances to an XYZ file.

        Args:
            idx1_list (list[int]): List of indices to save.
        """
        # Do't use sys.stdout = sys.__stdout__ to replace it, test function need use it
        original_stdout = sys.stdout
        with open(self.__filename, "w") as f:
            sys.stdout = f
            self.method_print(idx1_list)
        sys.stdout = original_stdout

    def method_save_xyz_append(self, idx1_list: list) -> None:  # append to old xyz file
        """
        Append selected Geometry instances to an existing XYZ file.

        Args:
            idx1_list (list): List of indices to append.
        """

        # Do't use sys.stdout = sys.__stdout__ to replace it, test function need use it
        original_stdout = sys.stdout
        with open(self.__filename, "a") as f:
            sys.stdout = f
            self.method_print(idx1_list)
        sys.stdout = original_stdout

    def method_print(self, idx1_St: list[int]) -> None:
        """
        Print selected Geometry instances to stdout.

        Args:
            idx1_St (list[int]): List of indices to print.
        """

        if (idx1_St == []):
            idx0_St: list = [*range(len(self))]
        else:
            idx0_St: list = [x-1 for x in idx1_St]
        for key in idx0_St:
            print(self.Sts[key], end="")

    def method_comment_keep(self) -> None:
        """
        Update comment fields for all Geometry instances.
        """

        for St in self.Sts:
            if St.method_update_comment():
                St.method_update_comment()

    def method_comment_new(self) -> None:
        """
        Assign unique cluster indices to all Geometry instances.
        """

        for idx0, St in enumerate(self.Sts):
            if St.method_update_comment():
                St.method_update_comment()
            St.method_comment_new(idx0+1)

    def method_rewrite_comment(self) -> None:
        """
        Format and update comment fields for all Geometry instances.
        """

        for St in self.Sts:
            if St.method_rewrite_comment():
                St.method_rewrite_comment()

    def get_comment_energy(self) -> list[float]:
        """
        Retrieve energy values from all Geometry instances.

        Returns:
            list [float]: List of energy values in Eh units.
        """

        Res: list = []
        for St in self.Sts:
            if St.get_comment_energy():
                Res.append(St.get_comment_energy())
        return Res

    def method_ensoGenFlexible(self, args, thermo_list) -> npt.NDArray:
        """
        Generate thermodynamic data for all Geometry instances.

        Args:
            args: Command-line arguments.
            thermo_list: Thermodynamic data list.

        Returns:
            npt.NDArray: Structured array with thermodynamic properties.
        """

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
        if np_enso.dtype.names is not None:
            names_anmr = list(np_enso.dtype.names)
        return np_enso[names_anmr[:8]]
