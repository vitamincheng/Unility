#!/usr/bin/env python
#  Module anmr/ censo           [08.27.2024] vitamin.cheng@gmail.com
from __future__ import annotations
from typing import Self
import os
import sys
import numpy as np
import numpy.typing as npt
from icecream import ic
from pathlib import Path
from censo_ext.Tools.utility import IsExist, IsExist_return_bool
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

    def __init__(self, DirFileName: Path) -> None:
        """
        Initialize the Anmrrc object by reading and parsing a .anmrrc file.

        Args:
            DirFileName (Path): The path to the .anmrrc file to be parsed.
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
            >>> anmrrc_obj = Anmrrc(Path("data/anmrrc_file.anmrrc"))
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

        lines: list[str] = open(DirFileName, "r").readlines()

        # Dict of Atomic Numbers and Atomic label
        self.Nums_element: dict[int, str] = {
            1: 'H', 6: 'C', 9: 'F', 14: 'Si', 15: 'P'}
        # List of Acid Atom (like NH OH) No show in spectra and the list of the number is atomic number
        # NH OH is 7 and 8
        self.acid_atoms_NoShow: list[int] = []
        # the 4 line of .anmrrc
        # [atomic number] [calculated shielding valule of the reference molecule] [experimental shift] [active or not]
        self.anmrrc: list[list[float]] = []
        # the active element if self.anmrrc that last number is 1 not 0
        self.Active: list[str] = []
        self.third_line: str = lines[2].rstrip()

        for idx, x in enumerate(lines[0].split()):
            if x != "XH":
                from censo_ext.Tools.utility import function_is_int
                if function_is_int(x):
                    self.acid_atoms_NoShow.append(int(x))
                else:
                    raise ValueError(
                        " Your .anmrrc file about 'XH acid atoms' havesomething wrong !!!")
            else:
                break

        _2nd_line: list[str] = lines[1].split()
        self.mf: float = float(_2nd_line[4])             # mf      : nmr frequency           # nopep8
        self.lw: float = float(_2nd_line[6])             # lw      : lines of width          # nopep8
        self.Temp: float = float(_2nd_line[12])          # Temp    : Temperature (K)         # nopep8
        self.JCoups: bool = bool(_2nd_line[8])           # JCoups  : bool of JCoups ONOFF    # nopep8
        self.SParams: bool = bool(_2nd_line[10])         # Sparams : bool of SParams ONOFF   # nopep8

        # 3 lines of .anmrrc Parameters
        for idx, x in enumerate(lines):
            if idx >= 3:
                self.anmrrc.append([int(x.split()[0]), float(
                    x.split()[1]), float(x.split()[2]), int(x.split()[3])])

        # Set the active species based on anmrrc entries
        for idx, x in enumerate(self.anmrrc):
            if x[3] == 1:
                self.Active.append(self.Nums_element[int(x[0])])

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

    def get_idx1_acid_atoms_NoShow_RemoveH(self, DirFileName: Path) -> list[int]:
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
            list[int]: A list of hydrogen atom indices (1-based) that should be removed 
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

        from censo_ext.Tools.Parameter import ELEMENT_NAMES
        acid_atoms_NoShow: list[str] = [ELEMENT_NAMES[i]
                                        for i in self.acid_atoms_NoShow]
        from censo_ext.Tools.ml4nmr import read_mol_neighbors
        mol, neighbors = read_mol_neighbors(DirFileName)
        acid_atoms_NoShowRemove: list[int] = []
        for idx, x in enumerate(mol):
            for y in acid_atoms_NoShow:
                if x.symbol == y:  # type: ignore
                    acid_atoms_NoShowRemove.append(idx+1)
        NoShow_Remove_Group: npt.NDArray[np.int64] = np.array(
            [], dtype=np.int64)
        for idx, x in enumerate(acid_atoms_NoShowRemove):
            NoShow_Remove_Group = np.concatenate(
                (NoShow_Remove_Group, neighbors[x]), axis=None)
        idx1_H_atom: list[int] = [idx+1 for idx,
                                 i in enumerate(mol) if i.symbol == "H"]  # type: ignore # nopep8
        return [i for i in NoShow_Remove_Group if i in idx1_H_atom]

    def get_Reference_anmrrc(self) -> float:
        """Retrieve the reference shielding value from the .anmrrc file.

        This method searches through the internal anmrrc data structure to find
        the reference shielding value. The reference is identified by looking
        for entries where the fourth element (index 3) equals 1.

        Returns:
            float: The reference shielding value found in the .anmrrc file.

        Raises:
            ValueError: If no active species (where x[3] == 1) is found in the
                .anmrrc file. This indicates that there is no valid reference
                shielding value available for processing.
        """

        reference: float | None = None
        for x in self.anmrrc:
            if x[3] == 1:
                reference = x[1]

        if reference:
            return reference
        else:
            print("No reference in your .anmrrc file")
            ic()
            raise ValueError("No reference in your .anmrrc file")


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

    def __init__(self, Directory: Path = Path(".")) -> None:
        """
        Initialize the Anmr object.

        Args:
            Directory (Path, optional): The directory containing NMR data files. Defaults to Path(".").
        """
        self.__Directory: Path = Directory
        self.enso: npt.NDArray                              # anmr_enso
        self.anmrJ: npt.NDArray[np.float64]                 # JCoup of anmr.out generated from anmr # nopep8
        self.anmrS: list[list[float]] = []                  # Shielding of anmr.out generated from anmr # nopep8
        # directory of orcaSJ
        self.orcaSJ: list[OrcaSJ] = []
        self.avg_orcaSJ = OrcaSJ()

        # anmr_nucinfo
        # idx1 and numbers of Chemical Equivalent
        self.nChemEqvs: dict[int, int] = {}
        # idx1 and neighbors index of Chemical Equivalent
        self.NeighborChemEqvs: dict[int, list[int]] = {}
        # idx1 and numbers of Magnetic Equivalent
        self.nMagnetEqvs: dict[int, int] = {}
        # idx1 and neighbors index of Magnetic Equivalent
        self.NeighborMangetEqvs: dict[int, list[int]] = {}

    def get_Directory(self) -> Path:
        """
        Get the directory path.

        Returns:
            Path: The directory containing NMR data files.
        """
        return self.__Directory

    def get_Anmr_Active(self) -> list[str]:
        """Get the list of active species from the .anmrrc file.

        This method retrieves the active atomic species that are currently
        configured in the anmrrc file. These species are typically used in
        ANMR calculations and simulations.

        Returns:
            list[str]: A list of atomic symbols (as strings) representing the
                active species defined in the .anmrrc configuration file.
        """
        return self.__AnmrParams.Active

    def get_idx1_acid_atoms_NoShow_RemoveH(self, DirFileName=Path("crest_conformers.xyz")) -> list[int]:
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

        return self.__AnmrParams.get_idx1_acid_atoms_NoShow_RemoveH(DirFileName)

    def method_read_anmrrc(self, fileName=Path(".anmrrc")) -> None:
        """Read .anmrrc setting file from censo.
        Args:
            fileName (Path, optional): Name of the .anmrrc file. Defaults to Path(".anmrrc").

        Raises:
            FileNotFoundError: If the specified .anmrrc file does not exist.
            Exception: If there are issues parsing the .anmrrc file.

        Example:
            >>> reader = AnmrFile()
            >>> reader.method_read_anmrrc(Path(".anmrrc"))
            >>> # File is read and parsed into self.__AnmrParams

        Note:
            The method sets the internal parameter `self.__AnmrParams` to an
            `Anmrrc` object created from the parsed file.
        """

        DirFileName: Path = self.__Directory / fileName
        IsExist(DirFileName)
        self.__AnmrParams: Anmrrc = Anmrrc(DirFileName)

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
            ic()
            raise FileNotFoundError(
                " Need to read the anmr_nucinfo and anmr_enso enso ")
        else:
            # for Normal of weight of anmr_enso
            weight: npt.NDArray[np.float64] = self.enso['BW']
            switch: npt.NDArray[np.int64] = self.enso['ONOFF']
            idx1_CONF: npt.NDArray[np.int64] = self.enso['CONF'].astype(
                np.int64)

            if np.sum(switch) == 0:
                print(" anmr_enso: Table - ONOFF is Zero ")
                ic()
                raise ValueError(" anmr_enso: Table - ONOFF is Zero ")

            weight = weight*switch
            weight = weight / np.sum(weight)
            normal_idx_weight: dict[np.int64, np.float64] = dict(
                zip(np.atleast_1d(idx1_CONF), np.atleast_1d(weight)))

            Active_orcaSJ: list[int] = []
            for idx, x in enumerate(self.orcaSJ):
                if x.CONFSerialNums in idx1_CONF:
                    Active_orcaSJ.append(idx)

            # orcaSParams and orcaJCoups using weighting to calculate and
            # save to Average_orcaSJ
            self.avg_orcaSJ: OrcaSJ = OrcaSJ()
            self.avg_orcaSJ.idx1Atoms = self.orcaSJ[0].idx1Atoms

            for idx in self.orcaSJ[0].SParams.keys():
                self.avg_orcaSJ.SParams[idx] = 0.0

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                idy0: list[int] = list(map(int, x.SParams.keys()))
                y: list[float] = list(map(float, x.SParams.values()))
                for idz, z in zip(idy0, np.array(y) * normal_idx_weight[x.CONFSerialNums]):
                    self.avg_orcaSJ.SParams[idz] += z.item()

            for key, value in self.avg_orcaSJ.SParams.items():
                self.avg_orcaSJ.SParams[key] = value - \
                    self.__AnmrParams.get_Reference_anmrrc()

            nShapes: int = np.shape(self.orcaSJ[0].JCoups[0])[0]
            self.avg_orcaSJ.JCoups = np.zeros((nShapes, nShapes))

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                y = x.JCoups
                self.avg_orcaSJ.JCoups += np.array(y) * \
                    normal_idx_weight[x.CONFSerialNums]

            ic(self.avg_orcaSJ.SParams)
            ic(self.avg_orcaSJ.JCoups)

        print(" ===== Finished the Average of all folder orcaS.out and orcaJ.out =====")

    def method_update_equiv_orcaSJ(self) -> None:
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
            ic()
            raise FileNotFoundError(
                " Need to read the anmr_nucinfo and anmr_enso enso ")
        else:
            print(" ===== Update the equivalent of SParams and JCoups =====")

            Atoms: list[int] = []
            for idx, x in self.orcaSJ[0].idx1Atoms.items():
                Atoms.append(idx)
            AtomsKeep: npt.NDArray[np.int64] = np.array(Atoms)
            AtomsEqvKeep: npt.NDArray[np.int64] = AtomsKeep.copy()

            ic(AtomsEqvKeep)

            for idx, x in enumerate(AtomsEqvKeep):
                if (self.nMagnetEqvs[x]) != 1:
                    if x != min(self.NeighborMangetEqvs[x]):
                        AtomsEqvKeep[idx] = 0
            AtomsEqvKeep = AtomsEqvKeep[np.nonzero(AtomsEqvKeep)]

            # Calculation the average ppm of Equivalent Atom and Replace the old ppm
            for orcaSJ in self.orcaSJ:
                for idx, x in enumerate(AtomsEqvKeep):
                    if (self.nChemEqvs[x] != 1):
                        ppm: list[float] = []
                        for y in self.NeighborChemEqvs[x]:
                            ppm.append(orcaSJ.SParams[y])
                        average: float = sum(ppm)/len(ppm)
                        for y in self.NeighborChemEqvs[x]:
                            orcaSJ.SParams[y] = average

            # for Equivalent Atom  of orcaJCoups
            # Calculation the average JCoups of Equivalent JCoups and Replace the old JCoups
            for orcaSJ in self.orcaSJ:
                for idx, x in enumerate(AtomsEqvKeep):
                    if (self.nMagnetEqvs[x] != 1):
                        for y in AtomsKeep:
                            JCoups: list[float] = []
                            average: float = 0
                            for z in (self.NeighborMangetEqvs[x]):
                                JCoups.append(orcaSJ.JCoups[list(
                                    AtomsKeep).index(y)][list(AtomsKeep).index(z)])
                                average: float = sum(JCoups)/len(JCoups)

                            for z in (self.NeighborMangetEqvs[x]):
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    y)][list(AtomsKeep).index(z)] = average
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    z)][list(AtomsKeep).index(y)] = average

                    if (self.nMagnetEqvs[x] > 2):
                        for k in (self.NeighborMangetEqvs[x]):
                            for ll in (self.NeighborMangetEqvs[x]):
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    k)][list(AtomsKeep).index(ll)] = 0
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    ll)][list(AtomsKeep).index(k)] = 0
                for idx, x in enumerate(orcaSJ.JCoups):
                    orcaSJ.JCoups[idx][idx] = 0

            AtomsDelete: list[int] = list(
                map(int, set(AtomsKeep).difference(set(list(AtomsEqvKeep)))))
            AtomsDelete.sort()

            # Delete Equivalent Atoms of idx1Atoms
            for orcaSJ in self.orcaSJ:
                for x in AtomsDelete:
                    del orcaSJ.idx1Atoms[x]

            # Delete Equivalent Atoms of orcaSParams
            for orcaSJ in self.orcaSJ:
                for x in AtomsDelete:
                    del orcaSJ.SParams[x]

            # Delete Equivalent Atom of orcaJCoups
            AtomsDelete2idx0: dict[int, int] = {}
            for idx0, x in enumerate(AtomsKeep):
                if x in AtomsDelete:
                    AtomsDelete2idx0[int(x)] = idx0
            ic(AtomsDelete2idx0)
            list_AtomsDelete: list[int] = [
                x for x in AtomsDelete2idx0.values()]
            list_AtomsDelete.reverse()

            for orcaSJ in self.orcaSJ:
                for x in list_AtomsDelete:
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, x, 0)
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, x, 1)

            idx1_acid_atoms_NoShow_RemoveH: list[int] = self.__AnmrParams.get_idx1_acid_atoms_NoShow_RemoveH(
                self.__Directory / Path("crest_conformers.xyz"))

            # Delete orcaSJ SParams in acid_atoms_NoShow
            idx0_AtomsDelete: list[int] = []
            for orcaSJ in self.orcaSJ:
                for idy0, SParam in enumerate(orcaSJ.SParams.copy()):
                    if SParam in idx1_acid_atoms_NoShow_RemoveH:
                        if idy0 not in idx0_AtomsDelete:
                            idx0_AtomsDelete.append(idy0)
                        del orcaSJ.SParams[SParam]

            # Delete orcaSJ orcaJCoups in acid_atoms_NoShow
            idx0_AtomsDelete.sort
            idx0_AtomsDelete.reverse()

            for orcaSJ in self.orcaSJ:
                for x in idx0_AtomsDelete:
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, x, 0)
                    orcaSJ.JCoups = np.delete(orcaSJ.JCoups, x, 1)

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

        Dir: Path = self.__Directory
        print("Files and directories in ", Dir, " : ")

        dirNames: list[str] = next(os.walk(Dir))[1]
        np.set_printoptions(formatter={'float': '{:12.5f}'.format})

        idx = 0
        while (idx != len(dirNames)):
            if (dirNames[idx].find('CONF') == -1):
                del dirNames[idx]
            else:
                idx += 1
        print("Directories = ", dirNames)
        del idx

        for idx0 in range(len(dirNames)):
            file_orcaS: Path = Dir / Path(dirNames[idx0] + "/NMR/orcaS.out")  # nopep8
            file_orcaJ: Path = Dir / Path(dirNames[idx0] + "/NMR/orcaJ.out")  # nopep8
            if os.path.exists(file_orcaS) and os.path.exists(file_orcaJ):
                print(str(idx0)+"  :  "+str(file_orcaS))
                print(str(idx0)+"  :  "+str(file_orcaJ))

                iter: OrcaSJ = OrcaSJ()
                iter.CONFSerialNums = int(dirNames[idx0].replace('CONF', ''))
                if not iter.method_read_orcaS(file=file_orcaS):
                    print("Something wrong in your orcaS.out")
                if not iter.method_read_orcaJ(file=file_orcaJ):
                    print("Something wrong in your orcaJ.out")
                self.orcaSJ.append(iter)

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

        fileName_Av_orcaS: Path = self.__Directory / \
            Path("Average/NMR/orcaS.out")
        fileName_Av_orcaJ: Path = self.__Directory / \
            Path("Average/NMR/orcaJ.out")
        fileName_Av_orcaAtoms: Path = self.__Directory / \
            Path("Average/NMR/orcaA.out")

        if IsExist_return_bool(fileName_Av_orcaS) and IsExist_return_bool(fileName_Av_orcaAtoms) and IsExist_return_bool(fileName_Av_orcaJ):
            return True
        else:
            return False

    def method_load_avg_orcaSJ(self, args) -> bool:
        """
        Load average orcaSJ data from files.

        This method reads and loads NMR data from ORCA calculation output files,
        including atom indices, scalar coupling constants (SParams), and J-couplings (JCoups).
        The method handles both BOBYQA and non-BOBYQA versions of the ORCA output files.

        Args:
            args: Command-line arguments containing BOBYQA flag.
                - args.bobyqa (bool): Flag indicating whether to use BOBYQA version
                    of the ORCA output files.

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
            >>> args = argparse.Namespace(bobyqa=True)
            >>> success = self.method_load_avg_orcaSJ(args)
            >>> print(success)
            True

        Files accessed:
            - Average/NMR/orcaS-BOBYQA.out (if args.bobyqa is True)
            - Average/NMR/orcaS.out (if args.bobyqa is False)
            - Average/NMR/orcaJ.out
            - Average/NMR/orcaA.out
        """

        from censo_ext.Tools.utility import jsonKeys2int
        if self.get_avg_orcaSJ_Exist():
            if args.bobyqa:
                file_avg_orcaS: Path = self.__Directory / \
                    Path("Average/NMR/orcaS-BOBYQA.out")
            else:
                file_avg_orcaS: Path = self.__Directory / \
                    Path("Average/NMR/orcaS.out")
            file_avg_orcaJ: Path = self.__Directory / \
                Path("Average/NMR/orcaJ.out")
            file_avg_orcaAtoms: Path = self.__Directory / \
                Path("Average/NMR/orcaA.out")

            import json
            with open(file_avg_orcaAtoms) as f:
                self.avg_orcaSJ.idx1Atoms = json.loads(
                    f.read(), object_pairs_hook=jsonKeys2int)

            from censo_ext.Tools.utility import load_dict_orcaS
            self.avg_orcaSJ.SParams = load_dict_orcaS(
                file_avg_orcaS)
            self.avg_orcaSJ.JCoups = np.loadtxt(file_avg_orcaJ)
            return True
        else:
            return False

    def method_save_adjust_avg_orcaS(self) -> None:
        """
        Save adjusted average orcaS data.

        Note: This method is currently not implemented.
        """
        raise NotImplementedError("Under Construct")

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

        avg_orcaS: Path = self.__Directory / Path("Average/NMR/orcaS.out")      # nopep8
        avg_orcaJ: Path = self.__Directory / Path("Average/NMR/orcaJ.out")      # nopep8
        avg_orcaAtoms: Path = self.__Directory / Path("Average/NMR/orcaA.out")  # nopep8

        Path(self.__Directory / Path("Average/NMR")
             ).mkdir(parents=True, exist_ok=True)

        from censo_ext.Tools.utility import save_dict_orcaS
        save_dict_orcaS(avg_orcaS, self.avg_orcaSJ.SParams)
        import json
        with open(avg_orcaAtoms, 'w') as f:
            f.write(json.dumps(self.avg_orcaSJ.idx1Atoms))
        np.savetxt(avg_orcaJ, self.avg_orcaSJ.JCoups, fmt="%10.5f")

    def method_save_folder_orcaSJ(self) -> None:
        """
        Save orcaSJ data to individual folder files.

        Note: This method is currently not implemented.
        """
        raise NotImplementedError("Under Construct")

    def method_read_anmrSJ(self, fileName: Path = Path("anmr.out")) -> None:
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

        fileName = self.__Directory / Path(fileName)
        IsExist(fileName)

        start_idx1: int = 0
        end_idx1: int = 0
        DataJ: list[str] = []
        lines: list[str] = open(fileName, "r").readlines()

        import re
        firstLine: bool = False
        start_shielding_idx1: int = 0
        for idx0, line in enumerate(lines):
            if re.search(r"MATRIX PRINTED:", line):
                start_idx1 = idx0 + 1
            if re.search(r"\+\/\-", line) and not firstLine:
                start_shielding_idx1 = idx0
                firstLine = True
        if start_shielding_idx1 != 0:
            nNuclei: int = start_idx1 - start_shielding_idx1 - 2
        else:
            print(" Something wrong ")
            ic()
            raise ValueError(" Something wrong in your anmr.out file")
        del firstLine, start_shielding_idx1

        nLines = 0
        for idx0 in range(int(nNuclei/6)):
            nLines += nNuclei - idx0 * 6 + 3
        end_idx1 = start_idx1 + nLines + nNuclei % 6 + 3 - 1

        for idx0, line in enumerate(lines):
            if idx0 >= start_idx1 and idx0 <= end_idx1:
                DataJ.append(line.rstrip())

        for line in lines:
            if re.search(r"\+\/\-", line):
                tmp: list[float] = [int(i) for i in line.split()[0:3]]
                tmp.append(float(line.split()[3]))
                self.anmrS.append(tmp)

        for line in lines:
            if re.search(r"1H resonance frequency", line):
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
        for idx, x in enumerate(self.anmrS):
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
            for y in self.NeighborChemEqvs[idx1]:
                print(f" {y:d}", end="")
            print("")
        for idx1 in self.nChemEqvs.keys():
            print(f"   {idx1:9d}   {self.nMagnetEqvs[idx1]:9d}")
            for y in self.NeighborMangetEqvs[idx1]:
                print(f" {y:4d}", end="")
            print("")

    def method_read_nucinfo(self, file: Path = Path("anmr_nucinfo")) -> None:
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

        file = self.__Directory / Path(file)
        IsExist(file)

        lines: list[str] = open(file, "r").readlines()
        del lines[0]

        Chemlines: list[str] = lines[0:int(len(lines)/2)]
        for idx, x in enumerate(Chemlines):
            x: str = x.rstrip()
            if (idx % 2) == 0:
                self.nChemEqvs[int(x.split()[0])] = int(x.split()[1])
            else:
                int_tmp: list[int] = []
                for y in x.split():
                    int_tmp.append(int(y))
                self.NeighborChemEqvs[int(x.split()[0])] = int_tmp

        Magnetlines: list[str] = lines[int(len(lines)/2):len(lines)]
        for idx, x in enumerate(Magnetlines):
            x: str = x.rstrip()
            if (idx % 2) == 0:
                self.nMagnetEqvs[int(x.split()[0])] = int(x.split()[1])
            else:
                int_tmp: list[int] = []
                for y in x.split():
                    int_tmp.append(int(y))
                self.NeighborMangetEqvs[int(x.split()[0])] = int_tmp

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
        """

        if len(in_np.dtype) != 8:  # type:ignore
            print("something wrong in your anmr_enso file")
            ic()
            raise ValueError("something wrong in your anmr_enso file")

    def method_read_enso(self, file: Path = Path("anmr_enso")) -> None:
        """Read ENSO data from file.

        This method reads ENSO data from a specified file
        and stores it in the instance variable `self.enso`. The file is expected to contain
        8 columns of data with specific meanings.

        Args:
            file (Path, optional): Name of the ENSO file to read. Defaults to Path("anmr_enso").
                The file is expected to be located in the directory specified by 
                `self.__Directory`.

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

        file = self.__Directory / Path(file)
        IsExist(file)

        self.enso = np.genfromtxt(file, names=True)
        if len(self.enso.dtype) != 8:                                       # type:ignore
            print("something wrong in your anmr_enso file")
            ic()
            raise FileNotFoundError("something wrong in your anmr_enso file")

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
        for Enso in self.enso:
            print(f'{int(Enso[0]):<1d}      {int(Enso[1]):<4d} {int(Enso[2]):<4d} {Enso[3]:>6.4f} {Enso[4]: > 11.7f} {Enso[5]: > 10.7f} {Enso[6]: > 10.7f} {Enso[7]:>4.3f}')  # nopep8

    def method_save_enso(self, file: Path = Path("anmr_enso.new")) -> None:
        """Save ENSO data to a file.

        This method writes the ENSO data to a specified file.
        The output is formatted using the print method for ENSO data.

        Args:
            file (Path, optional): The path to the output file. If not provided, 
                defaults to "anmr_enso.new" in the current directory.

        Example:
            >>> anmr = AnmrFile()
            >>> anmr.method_save_enso(Path("output_enso.txt"))
            >>> anmr.method_save_enso()  # Uses default filename

        Note:
            The method temporarily redirects stdout to the file during execution
            and then restores the original stdout.
        """

        DirFileName: Path = self.__Directory / Path(file)
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
        """
        self.JCoups: npt.NDArray[np.float64]
        self.SParams: dict[int, float] = {}
        self.Anisotropy: dict[int, float] = {}
        self.CONFSerialNums: int
        self.idx1Atoms: dict[int, str] = {}

    def method_read_orcaJ(self, file: Path = Path("orcaJ.out")) -> bool:
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

        print(" method_read_orcaJ", file)
        IsExist(file)

        start_idx: int
        end_idx: int
        start_idx, end_idx = 0, 0
        Data_str: list[str] = []
        DataJ: list[list[str]] = []
        lines: list[str] = open(file, "r").readlines()
        import re
        nVersion: str = ""

        for line in lines:
            if re.search(r"Program Version", line):
                nVersion: str = line

        if int(nVersion.split()[2][0]) == 5:
            nLines = 0
            start_idx = 0
            for idx, line in enumerate(lines):
                if re.search(r"Number of nuclei for epr/nmr", line):
                    nNuclei = int(line.split()[-1])
                    import math
                    nLines = (math.ceil(nNuclei/6))*(nNuclei+1)

                if re.search(r"SUMMARY OF ISOTROPIC COUPLING CONSTANTS", line):
                    start_idx = idx + 2

            if nLines != 0 and start_idx != 0:
                end_idx = start_idx + nLines - 1

        elif int(nVersion.split()[2][0]) == 6:
            for idx, line in enumerate(lines):
                if re.search(r"Maximum memory used throughout the entire PROP", line):
                    end_idx = idx - 4
                if re.search(r"SUMMARY OF ISOTROPIC COUPLING CONSTANTS", line):
                    start_idx = idx + 2
        else:
            print("This program is not work with before orca 5.0 ")

        if start_idx == 0 or end_idx == 0:
            print(file, " the data of the file is some error ...")
            print("  Exit and Close the program !!!")
            ic()
            raise ValueError(" the data of the file is some error ...")

        for idx, line in enumerate(lines):
            if idx >= start_idx and idx <= end_idx:
                Data_str.append(line.rstrip())

        for x in Data_str:
            DataJ.append(x.split()[2:])

        nums = 1
        while (nums >= 1):
            if ((int((nums-1)/6)+1)*(nums+1) == len(DataJ)):
                break
            nums: int = nums+1
        nAtomDataJ: int = nums

        for idx, x in enumerate(DataJ):
            if (idx > nAtomDataJ):
                DataJ[idx % (nAtomDataJ+1)] = DataJ[idx % (nAtomDataJ+1)]+x

        del DataJ[0]
        del DataJ[nAtomDataJ:]
        self.JCoups = np.array(DataJ, dtype=np.float64)
        return True

    def method_read_orcaS(self, file: Path = Path("orcaS.out")) -> bool:
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
            >>> success = reader.method_read_orcaS(Path("orcaS.out"))
            >>> print(success)
            True

        Note:
            This method populates the following instance attributes:
            - self.idx1Atoms: Dictionary mapping atom indices to atom symbols
            - self.SParams: Dictionary mapping atom indices to shielding constants
            - self.Anisotropy: Dictionary mapping atom indices to anisotropy values
        """

        print(" method_read_orcaS", file)
        IsExist(file)

        start_idx: int
        end_idx: int
        start_idx, end_idx = 0, 0
        DataS: list[str] = []
        lines: list[str] = open(file, "r").readlines()
        import re
        nVersion: str = ""
        nNuclei: int = 0

        for line in lines:
            if re.search(r"Program Version", line):
                nVersion: str = line

        if int(nVersion.split()[2][0]) == 5:
            for idx0, line in enumerate(lines):
                if re.search(r"Number of nuclei for epr/nmr", line):
                    nNuclei = int(line.split()[-1])
                if re.search(r"CHEMICAL SHIELDING SUMMARY", line):
                    start_idx = idx0 + 6
            end_idx = start_idx + nNuclei - 1
            if end_idx == 0 or start_idx == 0 or nNuclei == 0:
                return False
        elif int(nVersion.split()[2][0]) == 6:
            for idx0, line in enumerate(lines):
                if re.search(r"Maximum memory used throughout the entire PROP", line):
                    end_idx = idx0 - 5
                if re.search(r"CHEMICAL SHIELDING SUMMARY", line):
                    start_idx = idx0 + 6
        else:
            print(" This program is not work with before orca 5.0 ")
            ic()
            raise ValueError(" This program is not work with before orca 5.0 ")

        for idx0, line in enumerate(lines):
            if idx0 >= start_idx and idx0 <= end_idx:
                DataS.append(line.rstrip())
        self.idx1Atoms, self.Anisotropy, self.SParams = {}, {}, {}
        for x in DataS:
            self.idx1Atoms[int(x.split()[0])+1] = str(x.split()[1])
            self.SParams[int(x.split()[0])+1] = float(x.split()[2])
            self.Anisotropy[int(x.split()[0])+1] = float(x.split()[3])
        return True

    def method_save_orcaS(self) -> list:
        """
        Save orcaS data to list format.

        Note: This method is currently not implemented.
        """
        raise NotImplementedError("Under Construct")

    def method_save_orcaJ(self) -> list:
        """
        Save orcaJ data to list format.

        Note: This method is currently not implemented.
        """
        raise NotImplementedError("Under Construct")

    def method_print_orcaS(self) -> None:
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

        if len(self.idx1Atoms) == len(self.SParams):
            print("Nucleus  Element   Anisotropy")
            for idx, Atom in self.idx1Atoms.items():
                print(f'{idx:>5d}', f'{Atom:>8s}', end="")
                print(f'{self.SParams[idx]:>15.3f}')
        else:
            print("your orcaJ and orcaS is not fit each other")
            print("  Exit and Close the program !!!")
            ic()
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

        for idx0 in range(self.JCoups[0].size):
            for idy0 in range(self.JCoups[0].size):
                print(f'{(self.JCoups[idx0][idy0]):>8.3f}', end="")
            print("")
