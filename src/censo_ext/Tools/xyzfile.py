#!/usr/bin/env python

#  Modification             [08.26.2024] vitamin.cheng@gmail.com
#

from __future__ import annotations
from pathlib import Path
from typing import Self
import sys
import numpy as np
import numpy.typing as npt
import argparse
# from icecream import ic
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
        self.mass: npt.NDArray[np.float64]
        self.extras: list[list[str]] = extras
        self.com: npt.NDArray[np.float64]
        self.inertia: npt.NDArray[np.float64]

    def __add__(self, other: Self) -> Geometry:
        """ 
        Concatenate two Geometry objects by combining atoms and coordinates.

        Returns:
            Geometry: New Geometry instance with combined data.
        """

        import copy
        reGeometry: Geometry = copy.deepcopy(self)
        for key in other.names.keys():
            reGeometry.names[len(self.names) + key] = other.names[key]
        # coordinates of every atom
        reGeometry.coord = self.coord + other.coord
        reGeometry.extras = self.extras + other.extras
        reGeometry.nAtoms = len(self.names) + len(other.names)
        return reGeometry

    def __repr__(self) -> str:
        """ 
        Generate a string representation of the Geometry object in XYZ format.

        Returns:
            str: Formatted XYZ string with atom names, coordinates, and extras.
        """

        reStr: str = ""
        reStr += f"{self.nAtoms}\n{self.comment}\n"
        for idx0 in range(self.nAtoms):
            extra: str = "   ".join(
                self.extras[idx0]) if self.extras[idx0] != [] else ""
            reStr += f'{self.names[idx0+1]:>3s}    {self.coord[idx0][0]: 14.10f}'  # nopep8
            reStr += f'    {self.coord[idx0][1]: 14.10f}    {self.coord[idx0][2]: 14.10f}    {extra}\n'  # nopep8
        return reStr

    def get_comment_energy(self) -> float:
        """
        Retrieve the energy value from the comment field.

        Returns:
            float: Energy value in Eh units.
        """

        return self.comment_energy

    def method_translate_xyz(self, delta: npt.NDArray[np.float64]) -> Geometry:
        """
        Translate all atoms by a given vector.

        Args:
            delta (npt.NDArray): Translation vector (x, y, z).

        Returns:
            Geometry: Translated Geometry instance.
        """

        reGeometry: Geometry = copy.deepcopy(self)
        for x in reGeometry.coord:
            x += delta
        return reGeometry

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

    def method_idx_molecules_xyz(self, fileName: Path | str) -> list[list[int]]:
        """
        Identify molecular clusters from a topology file.

        Args:
            fileName (Path): Path to topology file.

        Returns:
            list[list[int]]: List of molecule atom index lists.
        """

        from censo_ext.Tools.topo import Topo
        import censo_ext.Tools.ml4nmr as ml4nmr
        fileName = Path(fileName)
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

    def method_computeCOM(self) -> npt.NDArray[np.float64]:
        """Calculate the center of mass (COM) of the molecule.

        This method computes the center of mass coordinates using the atomic masses
        and coordinates of all atoms in the molecule. If the center of mass has
        already been calculated and stored in self.com, it will be returned directly.
        Otherwise, it will compute the COM using the formula:
        COM = Σ(mass_i * coord_i) / Σ(mass_i)

        The method first updates the masses using method_update_masses() to ensure
        accurate calculations.

        Returns:
            npt.NDArray: An array of three floats representing the x, y, z coordinates
                of the center of mass in the same units as the atomic coordinates.
        """

        self.method_update_masses()
        if hasattr(self, 'com') and self.com.size == 3:
            return self.com
        else:
            self.com = np.dot(self.mass, self.coord) / np.sum(self.mass)
            return self.com

    def method_computeInertia(self) -> npt.NDArray[np.float64]:
        """Calculate the moment of inertia tensor.

        The moment of inertia tensor is computed using the standard formula:
        I_ij = -∑(m_k * (r_k,i * r_k,j - δ_ij * |r_k|^2))

        where:
        - m_k is the mass of atom k
        - r_k is the position vector of atom k relative to the center of mass
        - δ_ij is the Kronecker delta

        This implementation uses numpy's einsum for efficient computation.

        Returns:
            npt.NDArray: 3x3 inertia tensor representing the moment of inertia
                of the molecular system about the center of mass.
                The tensor is symmetric and represents the resistance to rotational motion
                around different axes.
        """

        data = self.coord - self.method_computeCOM()

        self.inertia = -np.einsum("ax,a,ay->xy", data, self.mass, data)
        return self.inertia

    def method_update_masses(self) -> None:
        """
        Update atomic masses based on element names.

        This method iterates through the atomic names stored in `self.names` and 
        retrieves the corresponding atomic weights from the `NAMES_WEIGHTS` 
        parameter dictionary. The weights are then appended to the `self.mass` array.

        Note:
            The element names in `self.names` should be in lowercase for proper 
            matching with the keys in `NAMES_WEIGHTS`.

        Example:
            >>> xyz = XYZFile()
            >>> xyz.names = {'atom1': 'C', 'atom2': 'H'}
            >>> xyz.method_update_masses()
            >>> print(xyz.mass)
            [12.011 1.008]

        Raises:
            KeyError: If an element name in `self.names` is not found in 
                      `NAMES_WEIGHTS`.
        """

        self.mass = np.array([])
        for value in self.names.values():
            from censo_ext.Tools.Parameter import NAMES_WEIGHTS
            self.mass = np.append(
                self.mass, NAMES_WEIGHTS[value.lower()], axis=None)

    def method_rewrite_comment(self) -> None:
        """
        Format and update the comment field with energy and cluster information.

        This method updates the object's comment field to include the energy value
        and cluster count in a standardized format. The energy is displayed with
        10 decimal places in Hartrees (Eh), and the cluster information is formatted
        as "#Cluster: {count}".

        The formatted comment follows this structure:
            " Energy =           {energy_value} Eh           #Cluster:     {cluster_count}"

        Example:
            If comment_energy = -76.45321 and comment_nClusters = 3,
            the result will be:
            " Energy =           -76.4532100000 Eh           #Cluster:     3"

        Note:
            This method modifies the object's comment field in-place.
        """

        self.comment = " Energy = "+" "*7 + \
            f"{self.comment_energy:.10f} Eh"+" "*8 + \
            f"#Cluster:     {self.comment_nClusters}"

    def method_update_comment(self) -> None:
        """
        Parse and extract energy and cluster information from the comment field.

        This method processes the comment line of an XYZ file to extract energy and 
        cluster information. It handles various formats of comment strings that may 
        contain energy values and cluster counts, and updates the object's 
        comment_energy and comment_nClusters attributes accordingly.

        The method expects comment strings in formats like:
        - "Energies= -76.432 #Cluster: 3"
        - "Energy= -76.432 #Cluster: 3"
        - "-76.432 #Cluster: 3"
        - "Eh -76.432 #Cluster: 3"

        If no valid energy or cluster information is found, default values of 0 
        are set for both attributes.

        Note:
            This method automatically calls method_rewrite_comment() after processing.

        Raises:
            ValueError: If the comment line cannot be parsed and contains invalid data.
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
            else:
                self.comment_energy, self.comment_nClusters = float(comments[0]), 0  # nopep8
        else:
            raise ValueError(
                f"{comments} Something wrong in your xyz file !!! ")

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

    def __init__(self, fileName: Path | str = Path("")) -> None:
        """
        Initialize a GeometryXYZs object to manage multiple Geometry instances.

        Args:
            fileName (Path): Path to XYZ file.
        """

        self.__filename: Path = Path(fileName)
        self.Sts: list[Geometry] = list()

    def __len__(self) -> int:
        return int(len(self.Sts))

    def set_filename(self, fileName: Path | str) -> None:
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

        xyzFile: GeometryXYZs = GeometryXYZs()
        if len(self) == 1:
            for x in np.linspace(0, 1, num=cut, endpoint=True):
                dSt = self.Sts[0].method_translate_xyz(delta=delta*x)
                xyzFile.Sts.append(copy.deepcopy(dSt))
            return xyzFile
        else:
            print(" More than 1 in your xyz file ")
            print(" All xyzs will to save your xyz file")
            for St in self.Sts:
                for x in np.linspace(0, 1, num=cut, endpoint=True):
                    dSt: Geometry = St.method_translate_xyz(delta=delta*x)
                    xyzFile.Sts.append(copy.deepcopy(dSt))
            return xyzFile

    def method_translate_xyzs(self, delta: npt.NDArray[np.float64]) -> GeometryXYZs:
        """
        Translate all Geometry instances in the collection.

        Args:
            delta (npt.NDArray[np.float64]): Translation vector.

        Returns:
            GeometryXYZs: Translated structures.
        """

        xyzFile: GeometryXYZs = GeometryXYZs()
        for x in self.Sts:
            x: Geometry = x.method_translate_xyz(delta=delta)
            xyzFile.Sts.append(copy.deepcopy(x))
        return xyzFile

    def __add__(self, Var: Self) -> GeometryXYZs:
        """
        Concatenate two GeometryXYZs objects.

        Args:
            Var (GeometryXYZs): Other GeometryXYZs instance.

        Returns:
            GeometryXYZs: Concatenated structures.
        """

        xyzFile: GeometryXYZs
        if len(Var) == 1:
            xyzFile = copy.deepcopy(self)
            xyzFile.Sts.append(Var.Sts[0])

            reFile: GeometryXYZs = GeometryXYZs()
            for idx in range(0, len(self), 1):
                reFile.Sts.append(
                    xyzFile.Sts[idx] + xyzFile.Sts[-1])
            return reFile
        else:
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
        for idx, x in enumerate(list_idx):
            if self.Sts[idx1+idx].method_molecules_separation_xyz(x):
                pass
            else:
                raise ValueError(
                    " Something wrong in your Molecule Separation xyz file")

        path: Path = Path("Separation")
        if not path.exists():
            path.mkdir()

        for plus_idx1, x in enumerate(list_idx, 1):
            self.set_filename(path / Path(str(plus_idx1)+".xyz"))
            self.method_save_xyz([idx1+plus_idx1])
        return True

    def method_Sts_extend(self, Sts_: list[Geometry]) -> None:
        raise NotImplementedError("Under Construct")
        # self.Sts.extend(Sts_)

    def method_Sts_append(self, St: Geometry) -> None:
        raise NotImplementedError("Under Construct")
        # self.Sts.append(St)

    def method_read_xyz(self) -> None:
        """Read XYZ file and populate the GeometryXYZs collection.

        This method reads a molecular geometry file in XYZ format and populates
        the collection with Geometry objects. Each frame in the XYZ file is
        converted to a Geometry object containing atomic names, coordinates,
        comment line, and any additional data fields.

        The XYZ file format expected is:
        - First line: number of atoms (nAtoms)
        - Second line: comment line
        - Next nAtoms lines: atom name and 3D coordinates (x, y, z)
        - Additional optional columns after coordinates are treated as extra data

        After reading all frames, the method calls method_comment_keep() to
        process the comment information.

        Note:
            The file must exist and be readable. The method uses IsExists_DirFileName
            to validate the file path before processing.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            ValueError: If the file format is invalid or contains parsing errors.
            IOError: If there are issues reading the file.

        Example:
            Given an XYZ file with content:
            2
            Water molecule
            O   0.000000   0.000000   0.000000
            H   0.757000   0.586000   0.000000

            The method will create a Geometry object with:
            - names: {1: 'O', 2: 'H'}
            - coords: [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0]]
            - comment: "Water molecule"
        """

        from censo_ext.Tools.utility import IsExists_DirFileName
        IsExists_DirFileName(self.__filename)

        with open(self.__filename, "r") as f:
            line: str = f.readline()
            while line != "":
                nAtoms = int(line)
                comment: str = f.readline().rstrip()
                names: dict[int, str] = dict()
                coords: list[npt.NDArray[np.float64]] = list()
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
        """Save selected Geometry instances to an XYZ file.

        Args:
            idx1_list (list[int]): List of indices to save. Each index corresponds to a
                Geometry instance that will be written to the output file.

        Example:
            >>> xyz_tool = XYZFileTool()
            >>> xyz_tool.method_save_xyz([0, 1, 2])
            # Saves the first three Geometry instances to the file

        Note:
            This method temporarily redirects sys.stdout to the output file during
            the printing process. The original stdout is restored after writing.
        """

        # Do't use sys.stdout = sys.__stdout__ to replace it, test function need use it
        original_stdout = sys.stdout
        with open(self.__filename, "w") as f:
            sys.stdout = f
            self.method_print(idx1_list)
        sys.stdout = original_stdout

    def method_save_xyz_append(self, idx1_list: list) -> None:  # append to old xyz file
        """Append selected Geometry instances to an existing XYZ file.

        Args:
            idx1_list (list): List of indices to append. Each index corresponds to a 
                Geometry instance that will be written to the XYZ file.
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

        This method iterates through all Geometry instances stored in self.Sts and
        collects their energy values. Each Geometry instance's get_comment_energy()
        method is called, and if it returns a valid energy value (not None), that
        value is added to the result list.

        The returned energies are in atomic units (Eh - Hartree).

        Returns:
            list[float]: A list of energy values in Hartree (Eh) units. 
                         Returns an empty list if no valid energies are found.
        """

        reEnergy: list = []
        for St in self.Sts:
            if St.get_comment_energy():
                reEnergy.append(St.get_comment_energy())
        return reEnergy

    def method_ensoGenFlexible(self, args, thermo_list) -> npt.NDArray:
        """Generate thermodynamic data for all Geometry instances.

        This method calculates various thermodynamic properties for a set of geometries
        including energies, Gibbs free energies, partition functions, and population
        distributions based on Boltzmann statistics.

        Args:
            args: Command-line arguments containing temperature information.
            thermo_list: List of thermodynamic data values for each geometry.

        Returns:
            npt.NDArray: Structured array with the following fields:
                - ONOFF: Flag indicating active status (int64)
                - NMR: NMR index (int64)
                - CONF: Conformer index (int64)
                - BW: Boltzmann weight (float64)
                - Energy: Electronic energy (float64)
                - Gsolv: Solvation free energy (float64)
                - mRRHO: Molecular RRHO correction (float64)
                - gi: Degeneracy factor (float64)

        The method performs the following calculations:
        1. Calculates total Gibbs free energy (Energy + mRRHO)
        2. Determines the minimum Gibbs free energy
        3. Computes delta Gibbs free energies relative to minimum
        4. Calculates partition functions (Qi) using Boltzmann distribution
        5. Computes population weights based on partition functions
        6. Normalizes weights to obtain final Boltzmann weights (BW)
        """

        import numpy.lib.recfunctions as rfn
        from censo_ext.Tools.Parameter import Eh, FACTOR
        # dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'),
        #       ('Energy', '<f8'), ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8')]

        # Column 8 is Total Gibbs Free Energy (Eh) = Energy + mRRHO
        TEMP: float = args.temp
        reEnso: npt.NDArray = np.zeros((len(self.Sts),), dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'),
                                                                ('Energy', '<f8'), ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8')])
        reEnso['ONOFF'] = 1
        reEnso['gi'] = 1.000
        reEnso['Gsolv'] = 0.00000000
        reEnso['NMR'] = np.arange(1, len(self.Sts)+1)
        reEnso['CONF'] = np.arange(1, len(self.Sts)+1)
        reEnso['mRRHO'] = np.array(thermo_list)
        reEnso['Energy'] = np.array(
            [a.comment_energy for a in self.Sts], dtype=[('Energy', 'f8')])
        Total: npt.NDArray = np.array(
            reEnso['Energy']+reEnso['mRRHO'], dtype=[('Total', 'f8')])
        reEnso = rfn.merge_arrays((reEnso, Total), flatten=True)

        # Gibbs_min is lowest energy of Gibbs Free Energy
        Gibbs_min: np.float64 = reEnso['Total'].min()

        # Column 9 is delta Gibbs Free Energy (kcal/mol)
        Gibbs: npt.NDArray = np.array(
            (reEnso['Total']-Gibbs_min)*Eh, dtype=[('Gibbs', 'f8')])
        reEnso = rfn.merge_arrays((reEnso, Gibbs), flatten=True)

        # Column 1o is Qi (each CONFS)
        Qi: npt.NDArray = np.array(
            np.exp(-reEnso['Gibbs']/(TEMP*FACTOR)), dtype=[('Qi', 'f8')])
        reEnso = rfn.merge_arrays((reEnso, Qi), flatten=True)

        # Qall is sum of Qi
        Qall: np.float64 = np.sum(reEnso['Qi'])

        # Column 11 is percentage of each CONFS
        NEW_BW: npt.NDArray = np.array(
            reEnso['Qi']/Qall, dtype=[('NEW_BW', 'f8')])
        reEnso = rfn.merge_arrays((reEnso, NEW_BW), flatten=True)
        # copy BW and delete residual parameter
        reEnso['BW'] = reEnso['NEW_BW']
        names_anmr: list = list()
        if reEnso.dtype.names:
            names_anmr = list(reEnso.dtype.names)
        return reEnso[names_anmr[:8]]
