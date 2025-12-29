#!/usr/bin/env python
from __future__ import annotations
from typing import Self
import sys
import re
import numpy as np
import numpy.typing as npt
from icecream import ic
from pathlib import Path
from censo_ext.Tools.utility import IntpID, IsExist, IsExist_bool, AtomID
from censo_ext.Tools.xyzfile import GeometryXYZs
# from dataclasses import dataclass


class Anmrrc():
    """Singleton class for handling .anmrrc files used in NMR calculations.

    This class reads and parses .anmrrc files that contain NMR reference parameters
    including atomic numbers, shielding values, experimental shifts, and active species.
    """

    _instance = None

    def __new__(cls, *args, **kwargs) -> Self:
        """
        Singleton Pattern
        """
        if not cls._instance:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, DirFile: Path | str) -> None:
        """
        Initialize the Anmrrc object by reading and parsing a .anmrrc file.

        Args:
            DirFile (Path | str): The path to the .anmrrc file to be parsed.
                This file contains NMR reference data including atomic numbers,
                calculated shielding values, experimental shifts, and active status.

        Attributes:
            Nums_element (dict[int, str]): Mapping of atomic numbers to their
                corresponding chemical symbols.
            acid_atoms_NoShow (list[int]): List of atomic numbers representing
                acid atoms (like NH, OH) that should not be shown in spectra.
            anmrrc (list[list[float]]): List of NMR parameters for each atom,
                where each entry contains [atomic_number, calculated_shielding,
                experimental_shift, active_status].
            Active (list[str]): List of chemical symbols for atoms marked as active.
            third_line (str): The content of the third line from the .anmrrc file.

        Example:
            >>> anmrrc_obj = Anmrrc("data/anmrrc_file.anmrrc")
            >>> print(anmrrc_obj.Active)
            ['H', 'C', 'F', 'Si', 'P']

        Note:
            The .anmrrc file format is expected to have:
            - First line: Acid atom identifiers (with "XH" as delimiter)
            - Second line: NMR parameters including frequency, linewidth, temperature,
              and flags for J-couplings and spin parameters
            - Third line: Additional metadata
            - Remaining lines: Atom-specific NMR data entries
        """
        DirFile = Path(DirFile)
        lines: list[str] = open(DirFile, "r").readlines()
        from censo_ext.Tools.utility import IsExists_DirFileName
        Dir, File = IsExists_DirFileName(DirFile)
        self.__Dir: Path = Path(Dir)
        # Dict of Atomic Numbers and Atomic label
        self.Nums_element: dict[int, str] = {
            1: 'H', 6: 'C', 9: 'F', 14: 'Si', 15: 'P'}
        # List of Acid Atom (like NH OH) No show in spectra and the list of the number is atomic number
        # NH OH is 7 and 8
        self.acid_atoms_NoShow: list[int] = []
        # the 4 line of .anmrrc
        # [atomic number] [calculated shielding valule of the reference molecule] [experimental shift] [active or not]
        # the active element if self.anmrrc that last number is 1 not 0
        self.anmrrc: list[list[float]] = []
        # the Active element is str : "H" or "C"
        self.Active: list[str] = []
        self.third_line: str = lines[2].rstrip()

        match: re.Match[str] | None = re.search(
            "(.*) XH acid atoms", lines[0])

        if match is None:
            pass
        else:
            for line in match.group(1).split():  # type: ignore
                from censo_ext.Tools.utility import function_is_int
                if function_is_int(line):
                    self.acid_atoms_NoShow.append(int(line))
                else:
                    raise ValueError(
                        " Your .anmrrc file about 'XH acid atoms' haves something wrong !!!")

        match = re.search(
            "mf= (.*) lw= (.*) J= (.*) S= (.*) T= (.*)", lines[1])

        self.mf: float = float(match.group(1))              # mf      : nmr frequency           # nopep8    # type: ignore
        self.lw: float = float(match.group(2))              # lw      : lines of width          # nopep8    # type: ignore
        self.JCoups: bool = bool(match.group(3))            # JCoups  : bool of JCoups ONOFF    # nopep8    # type: ignore
        self.SParams: bool = bool(match.group(4))           # Sparams : bool of SParams ONOFF   # nopep8    # type: ignore
        self.Temp: float = float(match.group(5))            # Temp    : Temperature (K)         # nopep8    # type: ignore

        # 3 lines of .anmrrc Parameters
        for line in lines[3:]:
            self.anmrrc.append([int(line.split()[0]), float(
                line.split()[1]), float(line.split()[2]), int(line.split()[3])])

        # Set the active species based on anmrrc entries
        for line in self.anmrrc:
            if line[3] == 1:
                self.Active.append(self.Nums_element[int(line[0])])

        # read the .anmrrc_linear the linear regression data
        self.linear: tuple[float, float] = self.get_anmrrc_linear()

    def __repr__(self) -> str:
        """
        Return a string representation of the Anmrrc object.

        Returns:
            str: Formatted string showing all parameters from the .anmrrc file.
        """
        Res: str = ""
        for x in (self.acid_atoms_NoShow):
            Res += f'{x} '
        Res += 'XH acid atoms\n'
        Res += f'ENSO qm= ORCA mf= {self.mf} lw= {self.lw}  J='
        Res += " on" if self.JCoups else " off"
        Res += " S="
        Res += ' on' if self.SParams else ' off'
        Res += f' T= {self.Temp}\n'
        Res += f'{self.third_line}\n'

        for x in self.anmrrc:
            Res += f"{x[0]:<d}  {x[1]:<9.2f} {x[2]:<6.1f} {x[3]:>2d}\n"
        return Res

    def get_acid_atoms_NoShow_RemoveH(self, DirFile: Path | str) -> list[AtomID]:
        """
        Get the list of hydrogen atom indices to be removed based on acid atom settings.

        This function identifies hydrogen atoms that should be removed from the molecular 
        structure based on the acid atom configuration. It first determines which atoms 
        are considered "acid atoms" according to the instance's acid_atoms_NoShow setting, 
        then finds all hydrogen atoms connected to these acid atoms and returns their indices.

        Args:
            DirFileName (Path): The path to the XYZ file containing the molecular structure.
                This file should contain the atomic coordinates and symbols for the molecule.

        Returns:
            list[AtomID]: A list of hydrogen atom indices (1-based) that should be removed 
                from the molecular structure. These are hydrogen atoms connected to the 
                acid atoms specified in acid_atoms_NoShow.

        Example:
            >>> remover = MoleculeRemover()
            >>> remover.acid_atoms_NoShow = [6, 8]  # Carbon and Oxygen
            >>> indices = remover.get_idx1_acid_atoms_NoShow_RemoveH(file_path)
            >>> print(indices)
            [5, 12, 18]  # Example hydrogen atom indices to remove

        Note:
            - The function uses 1-based indexing for atom indices as per standard molecular 
              file conventions
            - The returned indices are filtered to only include actual hydrogen atoms 
              present in the molecule
            - This function relies on the read_mol_neighbors utility to determine atomic 
              connectivity
        """
        DirFile = Path(DirFile)
        from censo_ext.Tools.Parameter import ELEMENT_NAMES
        # Element Names is str as "N" or "O" from 7 or 8
        acid_atoms_NoShow: list[str] = [ELEMENT_NAMES[i]
                                        for i in self.acid_atoms_NoShow]
        from censo_ext.Tools.ml4nmr import read_mol_neighbors
        xyzFile: GeometryXYZs = GeometryXYZs(DirFile)
        xyzFile.method_read_xyz()
        mol, neighbors = read_mol_neighbors(xyzFile=xyzFile)

        # find the atomID of molecule if is "N" or "O"
        acid_atoms_NoShowRemove: list[AtomID] = []
        for idx1, x in enumerate(mol, 1):
            if x.symbol in acid_atoms_NoShow:  # type: ignore
                acid_atoms_NoShowRemove.append(AtomID(idx1))

        NoShow_Remove: npt.NDArray[np.int64] = np.array(
            [], dtype=np.int64)
        for x in acid_atoms_NoShowRemove:
            NoShow_Remove = np.concatenate(
                (NoShow_Remove, neighbors[x]), axis=None)
        Has_H_atom: list[AtomID] = [idx1 for idx1,
                                 i in enumerate(mol, 1) if i.symbol == "H"]  # type: ignore # nopep8
        return [x for x in NoShow_Remove if x in Has_H_atom]

    def get_anmrrc_linear(self) -> tuple[float, float]:
        """Extract linear parameters from .anmrrc and .anmrrc_linear files.

        This method reads linear transformation parameters from the .anmrrc file and
        optionally from .anmrrc_linear file to determine the appropriate scaling
        and offset values for chemical shift calculations.

        The method searches for a reference point in the .anmrrc file (where column 3 equals 1)
        and uses it to determine whether to read from .anmrrc_linear or return a specific
        offset value. If a reference is found, it returns (-1, reference_value). Otherwise,
        it attempts to read two values from .anmrrc_linear file.

        Args:
            self: The instance of the class containing this method.

        Returns:
            tuple[float, float]: A tuple containing (a, b) where:
                - a: scaling factor for linear transformation
                - b: offset value for linear transformation

        Raises:
            SystemExit: If no reference is found in .anmrrc file or if .anmrrc_linear
                file is missing or contains invalid data.
        """
        reference: float | None = None
        for x in self.anmrrc:
            if x[3] == 1:
                reference = x[1]

        if reference is None:
            print("  no reference in your .anmrrc file")
            print("  exit and close the program !!!")
            exit(0)
        elif reference == 0:
            try:
                np_inData = np.genfromtxt(
                    self.__Dir / Path(".anmrrc_linear"), comments="#", usecols=[1])
            except FileNotFoundError:
                print("  no linear parameter in your .anmrrc_linear file")
                print("  exit and close the program !!!")
                exit(0)
            if len(np_inData) == 2:
                return (np_inData[0], np_inData[1])
            else:
                print(
                    "  linear parameters of your .anmrrc_linear file have Error, see the top words in .anmrrc_linear file")
                print("  exit and close the program !!!")
                exit(0)
        else:
            return (-1, reference)


class Anmr():
    """Main class for handling NMR data processing and analysis.

    This class serves as the central hub for processing NMR data from various sources
    including .anmrrc files, orcaS.out/orcaJ.out files, and anmr output files.
    """

    _instance = None

    def __new__(cls, *args, **kwargs) -> Self:
        '''Singleton Pattern'''
        if not cls._instance:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, Dir: Path | str = Path("."), verbose: bool = False) -> None:
        """
        Initialize the Anmr object.

        Args:
            Directory (Path | str): The directory containing NMR data files. Defaults to Path(".").
            verbose (bool, optional): Enable verbose output. Defaults to False.
        """
        self.__Dir: Path = Path(Dir)
        self.__verbose: bool = verbose
        self.enso: npt.NDArray                              # anmr_enso
        self.anmrJ: npt.NDArray[np.float64]                 # JCoup of anmr.out generated from anmr # nopep8
        self.anmrS: list[list[float]] = []                  # Shielding of anmr.out generated from anmr # nopep8
        # directory of orcaSJ
        self.orcaSJ: list[OrcaSJ] = []
        self.avg_orcaSJ = OrcaSJ()
        self.nNums_orcaS: int = 0
        self.nNums_orcaJ: int = 0

        # anmr_nucinfo
        # idx1 and numbers of Chemical Equivalent
        self.nChemEqvs: dict[AtomID, int] = {}
        # idx1 and neighbors index of Chemical Equivalent
        self.NeighborChemEqvs: dict[AtomID, list[AtomID]] = {}
        # idx1 and numbers of Magnetic Equivalent
        self.nMagnetEqvs: dict[AtomID, int] = {}
        # idx1 and neighbors index of Magnetic Equivalent
        self.NeighborMangetEqvs: dict[AtomID, list[AtomID]] = {}

        # For the data of Average Directory
        self.avg_Data_AD: Average_Directory = Average_Directory(self.__Dir)

    def get_Dir(self) -> Path:
        """
        Get the directory path.

        Returns:
            Path: The directory containing NMR data files.
        """
        return self.__Dir

    def get_Anmrrc_Active(self) -> list[str]:
        """Get the list of active species from the .anmrrc file.

        This method retrieves the active atomic species that are currently
        configured in the anmrrc file. These species are typically used in
        ANMR calculations and simulations.

        Returns:
            list[str]: A list of atomic symbols (as strings) representing the
                active species defined in the .anmrrc configuration file.
        """
        return self.__AnmrParams.Active

    def get_Anmrrc_linear(self) -> tuple[float, float]:
        return self.__AnmrParams.linear

    def get_idx1_acid_atoms_NoShow_RemoveH(self, DirFile: Path | str = Path("crest_conformers.xyz")) -> list[AtomID]:
        """
        Get the list of hydrogen atom indices to be removed based on acid atom settings.

        This function identifies hydrogen atoms that should be removed from the molecular 
        structure based on the acid atom configuration. It first finds all acid atoms 
        specified by the acid_atoms_NoShow attribute, then determines which hydrogen 
        atoms are bonded to these acid atoms and returns their indices.

        Args:
            DirFileName (Path): The path to the XYZ file containing the molecular structure.
                This file should contain the atomic coordinates and bonding information.

        Returns:
            list[int]: A list of hydrogen atom indices (1-based) that should be removed 
                from the molecular structure. These are hydrogen atoms bonded to the 
                acid atoms specified in self.acid_atoms_NoShow.

        Example:
            >>> remover = AtomRemover()
            >>> remover.acid_atoms_NoShow = [6, 8]  # Carbon and Oxygen
            >>> h_indices = remover.get_idx1_acid_atoms_NoShow_RemoveH("molecule.xyz")
            >>> print(h_indices)
            [5, 12, 18]

        Note:
            The function uses the molecular neighbors information from read_mol_neighbors 
            to determine which hydrogen atoms are bonded to the specified acid atoms.
        """
        DirFile = Path(DirFile)
        return self.__AnmrParams.get_acid_atoms_NoShow_RemoveH(DirFile)

    def method_read_anmrrc(self, file: Path | str = Path(".anmrrc")) -> None:
        """Read .anmrrc setting file from censo.

        This method reads and parses a .anmrrc configuration file, storing the
        parsed parameters in the internal `self.__AnmrParams` attribute as an
        `Anmrrc` object.

        Args:
            file (Path | str): Name of the .anmrrc file. Defaults to Path(".anmrrc").
                The file path is relative to the directory specified by `self.__Dir`.

        Example:
            >>> reader = AnmrFile()
            >>> reader.method_read_anmrrc(".anmrrc")
            >>> # File is read and parsed into self.__AnmrParams

        Note:
            The method sets the internal parameter `self.__AnmrParams` to an
            `Anmrrc` object created from the parsed file.
            Raises FileNotFoundError if the specified file does not exist.
        """
        file = Path(file)
        DirFile: Path = self.__Dir / file
        IsExist(DirFile)
        self.__AnmrParams: Anmrrc = Anmrrc(DirFile)

    def method_print_anmrrc(self) -> None:
        """Print the contents of the .anmrrc file.

        This method outputs the contents of the internal AnmrParams attribute
        to the standard output stream. The output includes all parameters
        stored in the .anmrrc configuration file.

        Args:
            self: The instance of the class containing this method.

        Returns:
            None: This method does not return any value.
        """

        print(self.__AnmrParams, end="")

    def method_avg_orcaSJ(self) -> None:
        """Average of orcaJ.out and orcaS.out in /NMR/CONFXX folder.

        This method computes the weighted average of NMR parameters (shielding constants
        and J-couplings) from multiple ORCA calculation folders. The averaging is
        performed using weights derived from the 'BW' column in the enso table, 
        adjusted by the 'ONOFF' flag to only include active configurations.

        The method populates the `avg_orcaSJ` attribute with the computed averages.
        It also subtracts a reference value from the shielding parameters using
        `get_Reference_anmrrc()`.

        Raises:
            FileNotFoundError: If nucinfo or enso data is not loaded.
            ValueError: If all entries in the 'ONOFF' column are zero.

        Note:
            The method assumes that self.orcaSJ contains valid OrcaSJ objects
            with properly initialized SParams and JCoups attributes.
        """

        print(" ===== Average of all folder orcaS.out and orcaJ.out =====")
        if len(self.nChemEqvs) == 0 or len(self.nMagnetEqvs) == 0 or self.enso.size == 0:
            print("  Need to read the anmr_nucinfo and anmr_enso enso ")
            print("  Exit and Close the program !!!")
            exit(0)
        else:
            # for Normal of weight of anmr_enso
            weight: npt.NDArray[np.float64] = self.enso['BW']
            switch: npt.NDArray[np.int64] = self.enso['ONOFF']
            idx1_CONFs: npt.NDArray[np.int64] = self.enso['CONF'].astype(
                np.int64)

            if np.sum(switch) == 0:
                print("  anmr_enso: Table - ONOFF is Zero ")
                print("  Exit and Close the program !!!")
                exit(0)

            weight = weight*switch
            weight = weight / np.sum(weight)
            normal_idx1_weight: dict[np.int64, np.float64] = dict(
                zip(np.atleast_1d(idx1_CONFs), np.atleast_1d(weight)))

            Active_orcaSJ: list[IntpID] = []
            for intp0, x in enumerate(self.orcaSJ):
                if x.CONF in idx1_CONFs:
                    Active_orcaSJ.append(IntpID(intp0))

            # orcaSParams and orcaJCoups using weighting to calculate and
            # save to Average_orcaSJ
            self.avg_orcaSJ: OrcaSJ = OrcaSJ()
            import copy
            self.avg_orcaSJ.Element = copy.deepcopy(self.orcaSJ[0].Element)

            # inital condition, let the chemical shift of average of orcaS is set to 0.0
            self.avg_orcaSJ.SParams = copy.deepcopy(self.orcaSJ[0].SParams)
            for key, _ in self.avg_orcaSJ.SParams.items():
                self.avg_orcaSJ.SParams[key] = 0.0

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                atomID: list[AtomID] = list(map(AtomID, x.SParams.keys()))
                ppm: list[float] = list(map(float, x.SParams.values()))
                for zip_atomID, weight_ppm in zip(atomID, np.array(ppm) * normal_idx1_weight[x.CONF]):
                    self.avg_orcaSJ.SParams[zip_atomID] += weight_ppm.item()

            # inital condition, let the JCoups of average of orcaSJ is set to 0.0 for every cell
            nShapes: int = np.shape(self.orcaSJ[0].JCoups[0])[0]
            self.avg_orcaSJ.JCoups = np.zeros((nShapes, nShapes))

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                self.avg_orcaSJ.JCoups += np.array(x.JCoups) * \
                    normal_idx1_weight[x.CONF]

            print("        Conf    Percentage(%)")
            for key, value in normal_idx1_weight.items():
                print(f"{key:12d} {value*100:12.3f}")

            if self.__verbose:
                ic(self.avg_orcaSJ.SParams)
                ic(self.avg_orcaSJ.JCoups)

        print(" ===== Finished the Average of all folder orcaS.out and orcaJ.out =====")

    def method_filter_active_orcaSJ(self, Active: str) -> None:
        """Filter ORCA SJ data to keep only atoms of specified active element.

        This method filters the ORCA SJ data structures to retain only atoms that
        belong to the specified active element. It removes atoms of other elements
        from Element, SParams, and JCoups dictionaries/arrays while maintaining the
        integrity of the remaining data structure.

        Args:
            self: The instance of the class containing this method.
            Active (str): The active element symbol to keep in the data structures.

        Returns:
            None: This method modifies the instance attributes orcaSJ in place.

        Note:
            The method prints filtering progress information to the console and
            performs in-place deletion of filtered atoms from all orcaSJ entries.
        """
        # Active Element - only for one element, "H" or "C"
        delete_atomIDs: list[AtomID] = [
            key for key, value in self.orcaSJ[0].Element.items() if value != Active]
        atomIDs: npt.NDArray[np.int64] = np.sort(
            np.array(list(self.orcaSJ[0].Element.keys())))
        delete_sorted_intp: npt.NDArray[np.intp] = np.array(
            atomIDs).searchsorted(delete_atomIDs)
        if len(delete_atomIDs) != 0:
            print(" ===== Filter the Active Atom of SParams and JCoups =====")
            for _orcaSJ in self.orcaSJ:
                for x in delete_atomIDs[::-1]:
                    if x in _orcaSJ.Element:
                        del _orcaSJ.Element[x]
                    if x in _orcaSJ.SParams:
                        del _orcaSJ.SParams[x]
                for x in delete_sorted_intp[::-1]:
                    _orcaSJ.JCoups = np.delete(_orcaSJ.JCoups, x, 0)
                    _orcaSJ.JCoups = np.delete(_orcaSJ.JCoups, x, 1)
            print(" ===== Finished the Filter of Active Atom of SParams and JCoups =====")

    def method_update_equiv_orcaSJ(self, Ref_FileName: Path = Path("crest_conformers.xyz")) -> None:
        """
        Update equivalent atoms in SParams and JCoups according to nuclear information.

        This method handles the replacement of equivalent atoms in both shielding parameters
        and coupling constants to ensure proper averaging across equivalent atoms. It processes
        the nuclear information to identify equivalent atoms and updates the corresponding
        SParams and JCoups values accordingly.

        The method performs the following operations:
        1. Validates that nucinfo and enso data are available
        2. Identifies equivalent atoms from the nuclear information
        3. Calculates average values for equivalent shielding parameters
        4. Calculates average values for equivalent coupling constants
        5. Removes duplicate equivalent atoms from the data structures
        6. Cleans up acid atoms from the final results

        Args:
            self: The instance of the class containing orcaSJ, nucinfo, and other relevant data

        Raises:
            FileNotFoundError: If nucinfo or enso data is not available

        Note:
            This method modifies the orcaSJ objects in-place by updating SParams and JCoups
            and removing equivalent atoms from idx1Atoms, SParams, and JCoups.
        """

        print(" Replace the equivalent of Sparams and JCoups")

        if len(self.nChemEqvs) == 0 or len(self.nMagnetEqvs) == 0 or self.enso.size == 0:
            print("  Need to read the anmr_nucinfo and anmr_enso enso ")
            print("  Exit and Close the program !!!")
            exit(0)
        else:
            print(" ===== Update the equivalent of SParams and JCoups =====")

            AtomIDs: list[AtomID] = [
                x for x in self.orcaSJ[0].Element.keys()]
            AtomIDsKeep: npt.NDArray[np.int64] = np.array(AtomIDs)
            AtomIDsEqvKeep: npt.NDArray[np.int64] = AtomIDsKeep.copy()

            if self.__verbose:
                ic(AtomIDsEqvKeep)

            for idx0, ids in enumerate(AtomIDsEqvKeep):
                if (self.nMagnetEqvs[ids]) != 1:
                    if ids != min(self.NeighborMangetEqvs[ids]):
                        AtomIDsEqvKeep[idx0] = 0
            AtomIDsEqvKeep = AtomIDsEqvKeep[np.nonzero(AtomIDsEqvKeep)]

            # Calculation the average ppm of Equivalent Atom and Replace the old ppm
            for orcaSJ in self.orcaSJ:
                for ids in AtomIDsEqvKeep:
                    if self.nChemEqvs[ids] != 1:
                        ppm: list[float] = []
                        for y in self.NeighborChemEqvs[ids]:
                            ppm.append(orcaSJ.SParams[y])
                        average: float = sum(ppm)/len(ppm)
                        for y in self.NeighborChemEqvs[ids]:
                            orcaSJ.SParams[y] = average

            # for Equivalent Atom of orcaJCoups
            # Calculation the average JCoups of Equivalent JCoups and Replace the old JCoups
            for orcaSJ in self.orcaSJ:
                for ids in AtomIDsEqvKeep:
                    if self.nMagnetEqvs[ids] != 1:
                        for y in AtomIDsKeep:
                            JCoups: list[float] = []
                            average: float = 0
                            intp_y = list(AtomIDsKeep).index(y)

                            for z in self.NeighborMangetEqvs[ids]:
                                intp_z = list(AtomIDsKeep).index(z)
                                JCoups.append(orcaSJ.JCoups[intp_y][intp_z])
                                average: float = sum(JCoups)/len(JCoups)

                            for z in self.NeighborMangetEqvs[ids]:
                                # intp_y=list(AtomIDsKeep).index(y)
                                intp_z = list(AtomIDsKeep).index(z)
                                orcaSJ.JCoups[intp_y][intp_z] = average
                                orcaSJ.JCoups[intp_z][intp_y] = average

                    if (self.nMagnetEqvs[ids] > 2):
                        for k, m in np.nditer([self.NeighborMangetEqvs[ids], self.NeighborMangetEqvs[ids]]):
                            intp_k = list(AtomIDsKeep).index(k)
                            intp_m = list(AtomIDsKeep).index(m)
                            orcaSJ.JCoups[intp_k][intp_m] = 0
                            orcaSJ.JCoups[intp_m][intp_k] = 0

                for idx0, ids in enumerate(orcaSJ.JCoups):
                    orcaSJ.JCoups[idx0][idx0] = 0

            AtomIDsDelete: list[AtomID] = list(
                map(AtomID, set(AtomIDsKeep).difference(set(list(AtomIDsEqvKeep)))))
            AtomIDsDelete.sort()

            # Delete Equivalent Atoms of idx1Atoms and orcaSParams
            for orcaSJ in self.orcaSJ:
                for ids in AtomIDsDelete:
                    del orcaSJ.Element[ids]
                    del orcaSJ.SParams[ids]

            # Delete Equivalent Atom of orcaJCoups
            DeleteAtomIDs2intpID: dict[AtomID, IntpID] = {}
            for idx0, ids in enumerate(AtomIDsKeep):
                if ids in AtomIDsDelete:
                    DeleteAtomIDs2intpID[ids] = IntpID(idx0)

            if self.__verbose:
                ic(DeleteAtomIDs2intpID)

            IntpID_Delete: list[IntpID] = [
                x for x in DeleteAtomIDs2intpID.values()]
            IntpID_Delete.reverse()

            for orcaSJ in self.orcaSJ:
                for ids in IntpID_Delete:
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, ids, 0)
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, ids, 1)

            Acid_atoms_NoShow_RemoveH: list[AtomID] = self.__AnmrParams.get_acid_atoms_NoShow_RemoveH(
                self.__Dir / Ref_FileName)
            print("="*80)

            # Delete orcaSJ SParams in acid_atoms_NoShow
            IntpID_Delete: list[IntpID] = []
            for orcaSJ in self.orcaSJ:
                for idy0, atomID in enumerate(orcaSJ.SParams.copy().keys()):
                    if atomID in Acid_atoms_NoShow_RemoveH:
                        if IntpID(idy0) not in IntpID_Delete:
                            IntpID_Delete.append(IntpID(idy0))
                        del orcaSJ.SParams[atomID]
                        del orcaSJ.Element[atomID]

            # Delete orcaSJ orcaJCoups in acid_atoms_NoShow
            IntpID_Delete.sort
            IntpID_Delete.reverse()

            for orcaSJ in self.orcaSJ:
                for ids in IntpID_Delete:
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, ids, 0)
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, ids, 1)

            print(" ===== Finished the equivalent of SParams and JCoups ===== ")

    def method_read_folder_orcaSJ(self) -> None:
        """Read orcaS.out and orcaJ.out files from all CONFXX directories.

        This method scans through all directories in the working directory that match
        the CONFXX pattern and reads NMR data from orcaS.out and orcaJ.out files.
        For each CONFXX directory, it attempts to read both orcaS.out and orcaJ.out files
        located in the NMR subdirectory. If both files exist, it processes them using
        OrcaSJ class methods and appends the results to self.orcaSJ list.

        Args:
            self: The instance of the class containing this method

        Returns:
            None: This method modifies the instance's orcaSJ attribute in-place

        Example:
            >>> analyzer = NMRAnalyzer()
            >>> analyzer.method_read_folder_orcaSJ()
            # Reads all orcaS.out and orcaJ.out files from CONF* directories
            # and populates analyzer.orcaSJ with parsed data

        Note:
            - Only directories matching the CONFXX pattern are processed
            - Files must exist in the NMR subdirectory of each CONFXX directory
            - Parsing errors in individual files will print error messages but won't stop processing
        """

        Dir: Path = self.__Dir
        print(f"Files and directories in {Dir} : ")
        dirNames: list[str] = [x for x in Dir.walk()][0][1]
        np.set_printoptions(formatter={'float': '{:12.5f}'.format})

        idx = 0
        while (idx != len(dirNames)):
            if (dirNames[idx].find('CONF') == -1):
                del dirNames[idx]
            else:
                idx += 1

        print(f"Directories = {dirNames}")
        del idx

        if len(dirNames) == 0:
            print("  Your CONFXX is not Exist !!!")
            print("  Exit and Close the program !!!")
            exit(0)
        from tqdm import tqdm
        for idx1, name in enumerate(tqdm(dirNames), 1):
            _orcaS: Path = Dir / Path(name + "/NMR/orcaS.out")  # nopep8
            _orcaJ: Path = Dir / Path(name + "/NMR/orcaJ.out")  # nopep8
            if self.__verbose:
                print(f"{idx1}  :  {_orcaS}")
                print(f"{idx1}  :  {_orcaJ}")

            iter: OrcaSJ = OrcaSJ()
            iter.CONF = int(name.replace('CONF', ''))
            if not iter.method_read_orcaS(file=_orcaS):
                print(" Your orcaS.out is missing or broken")
            else:
                self.nNums_orcaS += 1
            from censo_ext.Tools.utility import IsExist_bool
            if IsExist_bool(_orcaJ):
                iter.method_read_orcaJ(file=_orcaJ)
                self.nNums_orcaJ += 1
            else:
                print(" Your orcaJ.out is missing or broken")
            self.orcaSJ.append(iter)
        if self.nNums_orcaS == len(dirNames) and self.nNums_orcaJ == len(dirNames):
            return
        else:
            print(
                f"  {self.nNums_orcaS=}\n  {self.nNums_orcaJ=}\n  {len(dirNames)=}")
            print("  Your NMR folder haven't the same numbers to folder numbers !!!")
            print("  Exit and Close the program !!!")
            exit(0)

    # def get_avg_orcaSJ_Exist(self, AD: Average_Directory) -> bool:
    def get_avg_orcaSJ_Exist(self) -> bool:
        """
        Check if average orcaSJ files exist in the Average/NMR directory.

        This method verifies the existence of three specific ORCA output files
        in the Average/NMR subdirectory:
        - orcaS.out (structure file)
        - orcaJ.out (J-coupling file) 
        - orcaA.out (atom information file)

        The method performs a logical AND operation on all three file existence checks,
        returning True only if ALL files exist, False otherwise.

        Returns:
            bool: True if all average orcaSJ files exist in the Average/NMR directory,
                  False if any of the files are missing.
        """
        return self.avg_Data_AD.Exist()

    def method_BOBYQA_load_avg_orcaSJ(self) -> bool:
        """
        Load average orcaSJ data from files.

        This method reads and loads NMR data from ORCA calculation output files,
        including atom indices, scalar coupling constants (SParams), and J-couplings (JCoups).
        The method handles both BOBYQA and non-BOBYQA versions of the ORCA output files.


        Returns:
            bool: True if successful, False otherwise.
                - Returns False if average orcaSJ data does not exist.
                - Returns True after successfully loading all data from files.

        Note:
            This method modifies the instance attributes:
            - self.avg_orcaSJ.idx1Atoms: Dictionary of atom indices
            - self.avg_orcaSJ.SParams: Scalar coupling parameters
            - self.avg_orcaSJ.JCoups: J-coupling values

        Example:
            >>> success = self.method_load_avg_orcaSJ(args)
            >>> print(success)
            True

        Files accessed:
            - Average/NMR/orcaS-BOBYQA.out 
            - Average/NMR/orcaS.out 
            - Average/NMR/orcaJ.out
            - Average/NMR/orcaA.out
        """
        AD: Average_Directory = self.avg_Data_AD
        Result: bool = AD.method_load_files()
        self.avg_orcaSJ.Element = AD.Element
        if isinstance(AD.ChemicalShifts, dict):
            self.avg_orcaSJ.ChemicalShifts = AD.ChemicalShifts
            a, b = self.get_Anmrrc_linear()
            self.avg_orcaSJ.SParams = {
                key: (value-b)/a for key, value in AD.ChemicalShifts.items()}
            self.avg_orcaSJ.ChemicalShifts = {}
        elif isinstance(AD.ChemicalShifts, np.ndarray):
            temp: dict[AtomID, float] = {AtomID(key): float(value)
                                         for key, value, _ in AD.ChemicalShifts}
            a, b = self.get_Anmrrc_linear()
            self.avg_orcaSJ.SParams = {
                key: (value-b)/a for key, value in temp.items()}
            self.avg_orcaSJ.ChemicalShifts = {}
        else:
            print("  The type of your SParams have something wrong !!!")
            print("  Exit and Close the program !!!")
            exit(1)
        self.avg_orcaSJ.JCoups = AD.JCoups
        return Result

    def method_save_adjust_avg_orcaS(self) -> None:
        """
        Save adjusted average orcaS data.

        Note: This method is currently not implemented.
        """
        raise NotImplementedError("Under Construct")

    def method_linear_orcaS(self, inSParams: dict[AtomID, float]) -> dict[AtomID, float]:
        """Apply linear transformation to ORCA S parameters.

        This method applies a linear transformation to the input S parameters using
        the linear coefficients stored in self.__AnmrParams.linear. The transformation
        is defined as: output = a * input + b, where a and b are the linear coefficients.

        Args:
            inSParams (dict[AtomID, float]): Dictionary mapping atom IDs to their
                corresponding S parameters that need to be transformed.

        Returns:
            dict[AtomID, float]: Dictionary containing the transformed S parameters
                with the same atom ID keys as the input dictionary.
        """
        a, b = self.__AnmrParams.linear
        outSParams: dict[AtomID, float] = {key: a*value + b for key,
                                           value in inSParams.items()}
        return outSParams

    def method_save_avg_orcaSJ(self) -> None:
        """
        Save average orcaSJ data to files in Average/NMR directory.

        This method saves the averaged shielding parameters, coupling constants, and
        atom indices to respective files for future use.

        Args:
            self: The instance of the class containing this method.

        Returns:
            None: This method does not return any value.

        Example:
            >>> instance.method_save_avg_orcaSJ()
            # Saves data to Average/NMR/orcaS.out, Average/NMR/orcaJ.out, and Average/NMR/orcaA.out

        Note:
            The method creates the necessary directory structure (Average/NMR) if it doesn't exist.
            The saved files contain:
            - orcaS.out: Averaged shielding parameters
            - orcaJ.out: Coupling constants
            - orcaA.out: Atom indices
        """
        self.avg_Data_AD.Element = self.avg_orcaSJ.Element

        self.avg_Data_AD.ChemicalShifts = self.method_linear_orcaS(
            self.avg_orcaSJ.SParams)

        self.avg_Data_AD.JCoups = self.avg_orcaSJ.JCoups
        self.avg_Data_AD.method_save_files()
        self.avg_Data_AD.SParams = self.avg_orcaSJ.SParams

    def method_save_folder_orcaSJ(self) -> None:
        """
        Save orcaSJ data to individual folder files.

        Note: This method is currently not implemented.
        """
        raise NotImplementedError("Under Construct")

    def method_read_anmrSJ(self, fileName: Path | str = Path("anmr.out")) -> None:
        """
        Read the file anmr.out from anmr program.

        This method parses the anmr.out file to extract shielding constants and coupling
        constants (J-couplings) from the ANMR program output. It processes the matrix
        data to construct a symmetric J-coupling matrix and extracts resonance frequencies.

        Args:
            fileName (Path, optional): Name of the anmr.out file. Defaults to Path("anmr.out").

        Raises:
            ValueError: If the anmr.out file format is invalid or cannot be parsed correctly.
            FileNotFoundError: If the specified file does not exist.

        Note:
            The method expects specific formatting in the anmr.out file, including
            "MATRIX PRINTED:" markers and "+/-" lines for shielding data.

        Example:
            >>> reader = AnmrFile()
            >>> reader.method_read_anmrSJ("anmr.out")
            >>> print(reader.anmrS)  # Prints shielding constants
            >>> print(reader.anmrJ)  # Prints J-coupling matrix

        Side effects:
            - Populates self.anmrS with shielding data
            - Populates self.anmrJ with J-coupling matrix
            - Sets self.frq with resonance frequency
        """
        fileName = self.__Dir / Path(fileName)
        IsExist(fileName)

        start_idx1: int = 0
        end_idx1: int = 0
        DataJ: list[str] = []
        lines: list[str] = open(fileName, "r").readlines()

        firstLine: bool = False
        start_shielding_idx1: int = 0
        for idx0, line in enumerate(lines):
            if r"MATRIX PRINTED:" in line:
                start_idx1 = idx0 + 1
            if r"+/-" in line and not firstLine:
                start_shielding_idx1 = idx0
                firstLine = True
        if start_shielding_idx1 != 0:
            nNuclei: int = start_idx1 - start_shielding_idx1 - 2
        else:
            raise ValueError(" Something wrong in your anmr.out file")
        del firstLine, start_shielding_idx1

        nLines = 0
        for idx0 in range(int(nNuclei/6)):
            nLines += nNuclei - idx0 * 6 + 3
        end_idx1 = start_idx1 + nLines + nNuclei % 6 + 3 - 1

        for x in range(start_idx1, end_idx1+1):
            DataJ.append(lines[x].rstrip())

        for line in lines:
            if r"+/-" in line:
                tmp: list[float] = [int(i) for i in line.split()[0:3]]
                tmp.append(float(line.split()[3]))
                self.anmrS.append(tmp)

        for line in lines:
            if r"1H resonance frequency" in line:
                self.frq = float(line.split()[6])

        ListDataJ: list[str] = DataJ
        nLinesDataJ: int = len(self.anmrS)
        DataJDeleteBlank: list[str] = []

        k: int = nLinesDataJ
        i: int = len(ListDataJ)

        # Delete the top serial numbers of JCoup
        while i >= 1:
            del ListDataJ[0:3]
            DataJDeleteBlank += ListDataJ[0:k]
            ListDataJ = ListDataJ[k:]
            i, k = i-6, k-6
        del k, i

        # Delete the right serial numbers of JCoup
        DataDelete_nAtoms: list[str] = []
        for line in DataJDeleteBlank:
            if (len(DataDelete_nAtoms) % nLinesDataJ == 0):
                for j in range(6*((int)(len(DataDelete_nAtoms)/nLinesDataJ))):
                    DataDelete_nAtoms.append("\n")
            DataDelete_nAtoms.append(line[6:])

        DataJ_triangle: list[str] = [""]*(nLinesDataJ)

        # Restruct to one all JCoup table
        for idx0, line in enumerate(DataDelete_nAtoms):
            DataJ_triangle[idx0 % nLinesDataJ] = DataJ_triangle[idx0 %
                                                                nLinesDataJ].rstrip("\n") + " " + line.rstrip("\n")
        self.anmrJ = np.zeros((nLinesDataJ, nLinesDataJ))

        # copy half to other half data on JCoup
        for idx0 in range(nLinesDataJ):
            for idy0 in range(idx0):
                self.anmrJ[idy0][idx0] = self.anmrJ[idx0][idy0] = float(
                    DataJ_triangle[idx0].split()[idy0])

    def method_print_anmrS(self) -> None:
        """Print the anmr shielding data.

        This method displays the nucleus indices, coordination numbers, and chemical shift values
        in a formatted table. The output includes:
        - Nucleus index (5-digit right-aligned)
        - Coordination number (9-digit right-aligned)
        - Nucleus type (9-digit right-aligned)
        - Chemical shift value in ppm (13-digit right-aligned with 3 decimal places)

        Example output format:
            #  in coord file  # nucs   delta(ppm)
               1           6           8        123.456
               2           7           9        234.567

        Args:
            self: The instance of the class containing the anmrS data.

        Returns:
            None: This method prints directly to stdout and does not return any value.
        """

        print("    #  in coord file  # nucs   delta(ppm)")
        for x in self.anmrS:
            print(f"{x[0]:>5d} {x[1]:>9d} {x[2]:>9d} {x[3]:>13.3f}")

    def method_print_anmrJ(self) -> None:
        """Print the anmr coupling constants matrix.

        Displays the JCoups matrix in a formatted manner.

        The output format is a grid where each element is right-aligned in a
        field of width 10 with 5 decimal places. The matrix is printed row by
        row, with each row on a separate line.

        Example:
            If self.anmrJ is a 3x3 matrix, the output will look like:
            1.23456  2.34567  3.45678
            4.56789  5.67890  6.78901
            7.89012  8.90123  9.01234

        Note:
            This method modifies the standard output by printing directly to
            stdout without returning any value.
        """

        for idx0 in range(self.anmrJ[0].size):
            for idy0 in range(self.anmrJ[0].size):
                print(f'{self.anmrJ[idx0][idy0]:>10.5f}', end="")
            print("")

    def method_print_nucinfo(self) -> None:
        """
        Print nuclear information data.

        This method outputs the nuclear information stored in `self.nucinfo` to the console.
        The output includes:
        - The total number of atoms (first line)
        - For each atom group, the atom index and nucleus count
        - The equivalent atom groups for each atom

        The format is as follows:
        - First line: Total number of atoms (right-aligned, 12 digits)
        - Subsequent lines:
            - Atom index (3 digits) and nucleus count (3 digits)
            - List of equivalent atom indices (space-separated)

        Example output:
            123456789012
               1    2
             1 2 3 4
               5    1
             5

        Note:
            The data is read from `self.nucinfo` which should contain atom information
            structured as a list of lists, where each inner list represents an atom group.
        """

        nAtoms: int = len(self.nChemEqvs.keys())
        print(f"{nAtoms:>12d}")

        for idx1 in self.nChemEqvs.keys():
            print(f"   {idx1:d}   {self.nChemEqvs[idx1]:d}")
            for idy1 in self.NeighborChemEqvs[idx1]:
                print(f" {idy1:d}", end="")
            print("")

        for idx1 in self.nChemEqvs.keys():
            print(f"   {idx1:9d}   {self.nMagnetEqvs[idx1]:9d}")
            for idy1 in self.NeighborMangetEqvs[idx1]:
                print(f" {idy1:4d}", end="")
            print("")

    def method_read_nucinfo(self, file: Path | str = Path("anmr_nucinfo")) -> None:
        """
        Read nuclear information from a specified file.

        This method reads nuclear information from a file located in the object's directory.
        The file is expected to contain atomic data in a specific format where each atom's
        information is represented by two lines: first line contains atomic number and mass,
        second line contains additional integer values.

        Args:
            file (Path, optional): Name of the nucinfo file to read. Defaults to Path("anmr_nucinfo").
                The file is expected to be located in the object's directory.

        Raises:
            FileNotFoundError: If the specified file does not exist.
            ValueError: If the file format is incorrect or contains invalid data.
            IOError: If there are issues reading the file.

        Example:
            >>> obj.method_read_nucinfo("my_nucinfo_file")
            # Reads nuclear information from 'my_nucinfo_file' in the object's directory

        Note:
            The method expects the file to have a specific structure where:
            - First line contains number of atoms (nAtoms)
            - Subsequent lines contain atom data in pairs of lines
            - Each pair represents one atom with its properties
        """

        file = self.__Dir / Path(file)
        IsExist(file)

        lines: list[str] = open(file, "r").readlines()
        del lines[0]

        Chemlines: list[str] = lines[0:int(len(lines)/2)]
        for idx0, x in enumerate(Chemlines):
            x: str = x.rstrip()
            if (idx0 % 2) == 0:
                self.nChemEqvs[AtomID(int(x.split()[0]))] = AtomID(
                    int(x.split()[1]))
            else:
                ATomID_tmp: list[AtomID] = []
                for y in x.split():
                    ATomID_tmp.append(AtomID(int(y)))
                self.NeighborChemEqvs[AtomID(int(x.split()[0]))] = ATomID_tmp

        Magnetlines: list[str] = lines[int(len(lines)/2):len(lines)]
        for idx0, x in enumerate(Magnetlines):
            x: str = x.rstrip()
            if (idx0 % 2) == 0:
                self.nMagnetEqvs[AtomID(int(x.split()[0]))] = AtomID(
                    int(x.split()[1]))
            else:
                ATomID_tmp: list[AtomID] = []
                for y in x.split():
                    ATomID_tmp.append(AtomID(int(y)))
                self.NeighborMangetEqvs[AtomID(int(x.split()[0]))] = ATomID_tmp

    def method_create_enso(self, in_np: npt.NDArray) -> None:
        """Validate the enso data structure from an input numpy array.

        Args:
            in_np (npt.NDArray): Input numpy array containing enso data.

        Raises:
            ValueError: If the input numpy array does not have the expected
                dtype structure with exactly 8 fields.

        Example:
            >>> import numpy as np
            >>> data = np.array([(1, 2, 3, 4, 5, 6, 7, 8)],
            ...                 dtype=[('field1', 'i4'), ('field2', 'i4'),
            ...                        ('field3', 'i4'), ('field4', 'i4'),
            ...                        ('field5', 'i4'), ('field6', 'i4'),
            ...                        ('field7', 'i4'), ('field8', 'i4')])
            >>> obj.method_create_enso(data)
            # dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'), ('Energy', '<f8'),
            # ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8')])
        """

        if len(in_np.dtype) != 8:  # type:ignore
            raise ValueError(" Something wrong in your anmr_enso file")
        else:
            self.enso = np.array(in_np, dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'),
                                               ('BW', '<f8'), ('Energy',
                                                               '<f8'), ('Gsolv', '<f8'),
                                               ('mRRHO', '<f8'), ('gi', '<f8')])

    def method_read_enso(self, file: Path | str = Path("anmr_enso")) -> None:
        """Read ENSO data from file.

        This method reads ENSO data from a specified file
        and stores it in the instance variable `self.enso`. The file is expected to contain
        8 columns of data with specific meanings.

        Args:
            file (Path, optional): Name of the ENSO file to read. Defaults to Path("anmr_enso").
                The file is expected to be located in the directory specified by 
                `self.__Dir`.

        Raises:
            FileNotFoundError: If the specified file does not exist or if the file does not
                contain exactly 8 columns of data as expected.
            Exception: If there are issues with reading the file or parsing the data.

        Note:
            The expected column structure is:
            - ONOFF: On/Off flag
            - NMR: Nuclear Magnetic Resonance values
            - CONF: Configuration identifier
            - BW: Bandwidth
            - Energy: Energy values
            - Gsolv: Solvation free energy
            - mRRHO: Modified Reduced Rigid Rotor Harmonic Oscillator
            - gi: General index or identifier

        Example:
            >>> reader = AnmrFile()
            >>> reader.method_read_enso("my_enso_file")
            >>> print(reader.enso)
        """

        file = self.__Dir / Path(file)
        IsExist(file)

        self.enso = np.genfromtxt(file, names=True)
        if len(self.enso.dtype) != 8:                                       # type:ignore
            print("  something wrong in your anmr_enso file")
            print("  Exit and Close the program !!!")
            exit(0)

    def method_print_enso(self) -> None:
        """Print ENSO data in a formatted table.

        This method displays the ENSO (Electronic Nuclear Spin Orbits) data table
        containing various molecular properties. The table includes the following columns:

        - ONOFF: On/off flag (1 for on, 0 for off)
        - NMR: Nuclear Magnetic Resonance value
        - CONF: Configuration number
        - BW: Bandwidth value
        - Energy: Molecular energy value
        - Gsolv: Solvation free energy
        - mRRHO: Modified Reduced Rigid Rotor Harmonic Oscillator
        - gi: Partition function value

        The output is formatted with specific column widths and decimal precision
        for optimal readability.

        Example:
            >>> anmr_file.method_print_enso()
            ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi     
            1      1234   5678   0.1234     12345.6789012    12345.6789012    12345.6789012    1.234

        Note:
            The data is printed to standard output and formatted according to
            the internal data structure of self.enso.
        """

        print("ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi     ")
        for Enso in np.atleast_1d(self.enso):
            print(f'{int(Enso[0]):<1d}      {int(Enso[1]):<4d} {int(Enso[2]):<4d} {Enso[3]:>6.4f} {Enso[4]: > 11.7f} {Enso[5]: > 10.7f} {Enso[6]: > 10.7f} {Enso[7]:>4.3f}')  # nopep8

    def method_save_enso(self, file: Path | str = Path("anmr_enso.new")) -> None:
        """Save ENSO data to a file.

        This method writes the ENSO data to a specified file.
        The output is formatted using the print method for ENSO data.

        Args:
            file (Path, optional): The path to the output file. If not provided, 
                defaults to "anmr_enso.new" in the current directory.

        Example:
            >>> anmr = AnmrFile()
            >>> anmr.method_save_enso("output_enso.txt")
            >>> anmr.method_save_enso()  # Uses default filename

        Note:
            The method temporarily redirects stdout to the file during execution
            and then restores the original stdout.
        """

        DirFileName: Path = self.__Dir / Path(file)
        with open(DirFileName, "w") as f:
            sys.stdout = f
            self.method_print_enso()
        sys.stdout = sys.__stdout__


class OrcaSJ():
    """Class for handling ORCA NMR output files (orcaS.out and orcaJ.out).

    This class manages data from ORCA's chemical shielding and coupling constant outputs,
    including parsing and storing the information in structured formats.
    """

    def __init__(self) -> None:
        """
        Initialize the OrcaSJ object.

        This class holds data from ORCA's chemical shielding and coupling constant outputs.

        Attributes:
            JCoups (npt.NDArray[np.float64]): Coupling constants data.
            SParams (dict[AtomID, float]): Shielding parameters.
            ChemicalShifts (dict[AtomID, float]): ChemicalShits
            Anisotropy (dict[AtomID, float]): Anisotropy values.
            CONF (int): Configuration serial numbers.
            Element (dict[AtomID, str]): Mapping of atom indices to atom names.
            linear (tuple[float,float]): linear regression
        """
        self.JCoups: npt.NDArray[np.float64]
        self.SParams: dict[AtomID, float] = {}
        self.ChemicalShifts: dict[AtomID, float] = {}
        self.Anisotropy: dict[AtomID, float] = {}
        self.CONF: int
        self.Element: dict[AtomID, str] = {}
        self.linear: tuple[float, float]

    def method_load_anmrrc_linear(self, linear: tuple[float, float]) -> None:
        self.linear = linear

    def method_read_orcaJ(self, file: Path | str = Path("orcaJ.out")) -> bool:
        """
        Read the ORCA J-coupling output file for the censo program.

        This method parses the orcaJ.out file to extract isotropic coupling constants
        and stores them in the object's JCoups attribute.

        Args:
            file (Path, optional): Path to the orcaJ.out file. Defaults to Path("orcaJ.out").

        Returns:
            bool: True if successful, False otherwise.

        Raises:
            ValueError: If the data in the file is corrupted or incompatible with the expected format.
            FileNotFoundError: If the specified file does not exist.

        Note:
            This method supports ORCA versions 5.0 and 6.0. For other versions,
            a warning message will be printed and the program will exit.

        Example:
            >>> reader = AnmrFile()
            >>> success = reader.method_read_orcaJ("path/to/orcaJ.out")
            >>> print(reader.JCoups)
        """

        # print(f" method_read_orcaJ {file}")
        if not IsExist_bool(file):
            nShapes = len(self.SParams)
            self.JCoups = np.zeros((nShapes, nShapes))
            return True

        start_idx: int
        end_idx: int
        start_idx, end_idx = 0, 0
        Data_str: list[str] = []
        DataJ: list[list[str]] = []
        lines: list[str] = open(file, "r").readlines()
        nVersion: str = ""

        for line in lines:
            if r"Program Version" in line:
                nVersion: str = line

        if int(nVersion.split()[2][0]) == 5:
            nLines = 0
            start_idx = 0
            for idx0, line in enumerate(lines):
                if r"Number of nuclei for epr/nmr" in line:
                    nNuclei = int(line.split()[-1])
                    nLines = int(np.ceil(nNuclei/6))*(nNuclei+1)

                if r"SUMMARY OF ISOTROPIC COUPLING CONSTANTS" in line:
                    start_idx = idx0 + 2

            if nLines != 0 and start_idx != 0:
                end_idx = start_idx + nLines - 1

        elif int(nVersion.split()[2][0]) == 6:
            for idx0, line in enumerate(lines):
                if r"Maximum memory used throughout the entire PROP" in line:
                    end_idx = idx0 - 4
                if r"SUMMARY OF ISOTROPIC COUPLING CONSTANTS" in line:
                    start_idx = idx0 + 2
        else:
            print("This program is not work with before orca 5.0 ")

        if start_idx == 0 or end_idx == 0:
            raise ValueError(
                f"{file}, the data of the file have some error ...")

        for x in range(start_idx, end_idx+1):
            Data_str.append(lines[x].rstrip())

        for x in Data_str:
            DataJ.append(x.split()[2:])

        nums = 1
        while (nums >= 1):
            if ((int((nums-1)/6)+1)*(nums+1) == len(DataJ)):
                break
            nums: int = nums+1
        nAtomDataJ: int = nums

        for idx0, x in enumerate(DataJ):
            if (idx0 > nAtomDataJ):
                DataJ[idx0 % (nAtomDataJ+1)] = DataJ[idx0 % (nAtomDataJ+1)]+x

        del DataJ[0]
        del DataJ[nAtomDataJ:]
        self.JCoups = np.array(DataJ, dtype=np.float64)
        return True

    def method_read_orcaS(self, file: Path | str = Path("orcaS.out")) -> bool:
        """
        Read the ORCA S (NMR shielding) output file for the censo program.

        This method parses the ORCA S output file to extract NMR shielding data,
        including isotropic shielding constants and anisotropy values for each nucleus.

        Args:
            file (Path, optional): Path to the orcaS.out file. Defaults to Path("orcaS.out").

        Returns:
            bool: True if successful, False otherwise.

        Raises:
            ValueError: If the ORCA version is not supported (must be version 5.0 or higher).

        Example:
            >>> reader = NMRFileReader()
            >>> success = reader.method_read_orcaS("orcaS.out")
            >>> print(success)
            True

        Note:
            This method populates the following instance attributes:
            - self.idx1Atoms: Dictionary mapping atom indices to atom symbols
            - self.SParams: Dictionary mapping atom indices to shielding constants
            - self.Anisotropy: Dictionary mapping atom indices to anisotropy values
        """
        IsExist(file)

        start_idx: int
        end_idx: int
        start_idx, end_idx = 0, 0
        DataS: list[str] = []
        lines: list[str] = open(file, "r").readlines()
        nVersion: str = ""
        nNuclei: int = 0

        for line in lines:
            if r"Program Version" in line:
                nVersion: str = line

        if int(nVersion.split()[2][0]) == 5:
            for idx0, line in enumerate(lines):
                if r"Number of nuclei for epr/nmr" in line:
                    nNuclei = int(line.split()[-1])
                if r"CHEMICAL SHIELDING SUMMARY" in line:
                    start_idx = idx0 + 6
            end_idx = start_idx + nNuclei - 1
            if end_idx == 0 or start_idx == 0 or nNuclei == 0:
                return False
        elif int(nVersion.split()[2][0]) == 6:
            for idx0, line in enumerate(lines):
                if r"Maximum memory used throughout the entire PROP" in line:
                    end_idx = idx0 - 5
                if r"CHEMICAL SHIELDING SUMMARY" in line:
                    start_idx = idx0 + 6
        else:
            print(" This program is not work with before orca 5.0 ")
            print("  Exit and Close the program !!!")
            exit(0)

        for x in range(start_idx, end_idx+1):
            DataS.append(lines[x].rstrip())

        self.Element, self.Anisotropy, self.SParams = {}, {}, {}
        for x in DataS:
            idx1: AtomID = AtomID(int(x.split()[0]) + 1)
            self.Element[idx1] = str(x.split()[1])
            self.SParams[idx1] = float(x.split()[2])
            self.Anisotropy[idx1] = float(x.split()[3])
        return True

    def method_setup_ChemicalShifts(self) -> None:
        """Setup chemical shifts using linear transformation.

        This method applies a linear transformation to the S parameters to calculate
        chemical shifts. It uses the linear coefficients stored in self.linear to
        perform the transformation: chemical_shift = a * s_param + b. The results are
        stored in self.ChemicalShifts dictionary.

        The method also prints diagnostic information about the linear regression
        coefficients to the console, including the equation format and coefficient values.

        Args:
            self: The instance of the class containing this method.

        Returns:
            None: This method modifies the instance attribute self.ChemicalShifts
                in place rather than returning a value.
        """
        a, b = self.linear
        print(" ===== Print the Linear Regression =====")
        print("  y = ax + b")
        print(f"  a = {a}     b = {b}")
        print("  [see .anmrrc and .anmrrc_linear]")
        print("")
        self.ChemicalShifts = {
            key: a*value + b for key, value in self.SParams.items()}

    def method_teardown_ChemicalShifts(self) -> None:
        self.ChemicalShifts = {}

    def method_print_av_orcaS(self) -> None:
        """Print ORCA-S data.

        This method displays the nucleus indices, element symbols, and chemical
        shielding values from the ORCA-S calculation results.

        The output format is:
        - Nucleus: atom index (5-digit right-aligned)
        - Element: element symbol (8-character right-aligned)
        - Anisotropy: chemical shielding value (15-character right-aligned, 3 decimal places)

        Raises:
            ValueError: If the number of atoms does not match the number of shielding parameters,
                       indicating a mismatch between ORCA-J and ORCA-S data.
            ic: If a mismatch is detected, the program will exit with an error message.

        Note:
            This method requires self.idx1Atoms and self.SParams to be properly initialized
            with matching lengths for correct operation.
        """
        if len(self.Element) == len(self.ChemicalShifts):
            print(" ===== Print the Chemical Shift of Atoms =====")
            print("    coord  Element     Anisotropy")
            for idx1, Element in self.Element.items():
                print(f'   {idx1:>5d}', f'{Element:>8s}', end="")
                print(f'{self.ChemicalShifts[idx1]:>15.3f}')
            print("")
        else:
            raise ValueError("your orcaJ and orcaS is not fit each other")

    def method_print_orcaJ(self) -> None:
        """Print orcaJ data.

        Displays the JCoups matrix in a formatted manner.

        The output format is a table where each row represents a nucleus and each
        column represents the coupling constant between that nucleus and all other
        nuclei. The values are printed with 3 decimal places in a fixed-width
        format.

        Example:
            If JCoups is a 3x3 matrix, the output will look like:
            10.500  2.300  0.800
             2.300 15.200  1.100
             0.800  1.100  8.700

        Note:
            This method modifies the standard output stream directly.
            The matrix is assumed to be square and symmetric (for coupling constants).
        """
        print(" ===== Print the Coupling Constant =====")
        for idx0 in range(self.JCoups[0].size):
            for idy0 in range(self.JCoups[0].size):
                print(f'{(self.JCoups[idx0][idy0]):>8.3f}', end="")
            print("")
        print("")


class Average_Directory(object):
    """Base class for handling average directory operations in ANMR calculations.

    This is the base class that provides common functionality for managing
    directories and files used in ANMR (Average Nuclear Magnetic Resonance)
    calculations. It establishes the directory structure and defines the
    necessary file paths and data attributes for ORCA calculation files.

    The class manages three primary file types:
    - ORCA S files (spin parameters)
    - ORCA J files (J coupling constants) 
    - ORCA A files (atomic indices)

    Attributes:
        _Dir: Path object representing the average NMR directory
        _orcaS: Path object for the ORCA S output file (not initialized in base class)
        _orcaJ: Path object for the ORCA J output file
        _orcaA: Path object for the ORCA A output file
        SParams: Dictionary or NumPy array containing spin parameters
        JCoups: NumPy array containing J coupling constants
        idx1Atoms: Dictionary containing atomic indices

    Args:
        Dir: Path object representing the base directory. Defaults to current directory (".")

    Example:
        >>> avg_dir = Average_Directory(Path("./analysis"))
        >>> print(avg_dir._Dir)
    """
    FileName_A = Path("orcaA.out")
    FileName_J = Path("orcaJ.out")
    FileName_Average = Path("Average/NMR")

    def __init__(self, Dir: Path = Path(".")) -> None:
        self._Dir: Path = Dir / self.FileName_Average
        self._orcaS: Path
        self._orcaJ: Path = self._Dir / self.FileName_J
        self._orcaA: Path = self._Dir / self.FileName_A
        self.SParams: dict[AtomID, float] | npt.NDArray
        self.ChemicalShifts: dict[AtomID, float] | npt.NDArray
        self.JCoups: npt.NDArray
        self.Element: dict[AtomID, str]

    def method_print_file(self) -> None:
        """Print the file paths of ORCA A, S, and J files."""
        print(self._orcaA)
        print(self._orcaS)
        print(self._orcaJ)

    def method_load_files(self) -> bool:
        """Load ORCA files for ANMR analysis.

        This method loads three types of files required for ANMR calculations:
        - ORCA A file (atomic indices)
        - ORCA S file (spin parameters)
        - ORCA J file (J coupling constants)

        The method performs the following operations:
        1. Checks if all three files exist
        2. Loads atomic indices from the ORCA A file using jsonKeys2int conversion
        3. Loads spin parameters from the ORCA S file, handling two possible formats:
           - Two-column format: integer index and float value
           - Three-column format: loaded with np.loadtxt
        4. Loads J coupling constants from the ORCA J file

        Args:
            self: The instance of the class containing this method

        Returns:
            bool: True if all files were successfully loaded, False otherwise

        Raises:
            SystemExit: If the SParams format is invalid (not 2 or 3 columns)
        """
        from censo_ext.Tools.utility import jsonKeys2int
        if self.Exist():   # all three files must be exists

            # load the orcaA file
            import json
            with open(self._orcaA) as f:
                self.Element = json.loads(
                    f.read(), object_pairs_hook=jsonKeys2int)

            lines: list = open(self._orcaS, "r").readlines()
            if len(lines[0].split()) == 2:
                Data: dict[AtomID, float] = {}
                for x in lines:
                    Data[AtomID(int(x.split()[0]))] = float(x.split()[1])
                self.ChemicalShifts = Data
            elif len(lines[0].split()) == 3:
                self.ChemicalShifts = np.loadtxt(self._orcaS)
            else:
                print("  The type of your SParams have something wrong !!!")
                print("  Exit and Close the program !!!")
                exit(1)

            # load the Joups file
            self.JCoups = np.loadtxt(self._orcaJ)
            return True
        else:
            return False

    def method_save_files(self) -> None:
        """Save ANMR files to disk.

        This method saves three types of files required for ANMR calculations:
        - ORCA A file (atomic indices in JSON format)
        - ORCA S file (spin parameters in formatted text)
        - ORCA J file (J coupling constants in formatted text)

        The method handles different data formats for spin parameters:
        - Dictionary format: saves as two-column text file
        - 2D NumPy array: saves with appropriate formatting based on column count
        - 3D NumPy array: saves with appropriate formatting based on column count

        Args:
            self: The instance of the class containing this method

        Returns:
            None: This method doesn't return anything, it performs file I/O operations

        Raises:
            SystemExit: If the SParams format is invalid or has incorrect dimensions
        """
        (self._Dir).mkdir(parents=True, exist_ok=True)
        if isinstance(self.ChemicalShifts, dict):
            with open(self._orcaS, 'w') as f:
                for key, value in self.ChemicalShifts.items():
                    f.write(f'{key:10d} {value:12.5f}\n')
        elif isinstance(self.ChemicalShifts, np.ndarray):
            if self.ChemicalShifts.shape[1] == 3:
                np.savetxt(self._orcaS, self.ChemicalShifts,
                           fmt="%10d   %10.5f %10d")
            elif self.ChemicalShifts.shape[1] == 2:
                np.savetxt(self._orcaS, self.ChemicalShifts,
                           fmt="%10d   %10.5f")
            else:
                print("  The type of your SParams.shape have someting wrong !!!!")
                print("  Exit and Close the program !!!")
                exit(1)

        else:
            print("  The type of your SParams have someting wrong !!!!")
            print("  Exit and Close the program !!!")
            exit(1)

        import json
        with open(self._orcaA, 'w') as f:
            f.write(json.dumps(self.Element))

        np.savetxt(self._orcaJ, self.JCoups, fmt="%10.5f")

    def Exist(self) -> bool:
        """Check if all required ORCA files exist.

        This method verifies that all three necessary files for ANMR calculations exist
        in the specified directory. The required files are:
        - ORCA A file (atomic indices)
        - ORCA S file (spin parameters)
        - ORCA J file (J coupling constants)

        Args:
            self: The instance of the class containing this method

        Returns:
            bool: True if all three files exist, False otherwise
        """

        if all(IsExist_bool(file_path) for file_path in (self._orcaS, self._orcaJ, self._orcaA)):
            return True
        else:
            return False


class AD_Normal(Average_Directory):
    """A class representing normal average directory for ANMR calculations.

    This class extends Average_Directory and is specifically designed for handling
    ORCA files in normal ANMR calculations. It manages the file paths for the three
    required ORCA files (S, A, and J) and provides methods for loading and saving
    these files.

    Attributes:
        _orcaS: Path object pointing to the ORCA S output file
        _orcaA: Path object pointing to the ORCA A output file  
        _orcaJ: Path object pointing to the ORCA J output file

    Args:
        Dir: Path object representing the directory containing the files.
             Defaults to current directory (".")

    Example:
        >>> ad = AD_Normal(Path("./my_analysis"))
        >>> ad.load_files()
    """
    FileName_S = Path("orcaS.out")

    def __init__(self, Dir: Path = Path(".")) -> None:
        super().__init__(Dir)
        self._orcaS = self._Dir / self.FileName_S


class AD_BOBYQA(Average_Directory):
    """A class representing BOBYQA average directory for ANMR calculations.

    This class extends Average_Directory and is specifically designed for handling
    ORCA files in ANMR calculations using the BOBYQA optimization method. It manages
    the file paths for the three required ORCA files (S, A, and J) with BOBYQA-specific
    naming conventions.

    Attributes:
        _orcaS: Path object pointing to the BOBYQA ORCA S output file
        _orcaA: Path object pointing to the BOBYQA ORCA A output file  
        _orcaJ: Path object pointing to the BOBYQA ORCA J output file

    Args:
        Dir: Path object representing the directory containing the files.
             Defaults to current directory (".")

    Example:
        >>> ad = AD_BOBYQA(Path("./bobyqa_analysis"))
        >>> ad.load_files()
    """

    FileName_S = Path("orcaS-BOBYQA.out")

    def __init__(self, Dir: Path = Path(".")) -> None:
        super().__init__(Dir)
        self._orcaS = self._Dir / self.FileName_S
