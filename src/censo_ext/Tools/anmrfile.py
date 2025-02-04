#!/usr/bin/env python3
#  Module anmr/ censo           [08.27.2024] vitamin.cheng@gmail.com
from __future__ import annotations
# from numpy import exceptions
from typing_extensions import Self
import os
import sys
import numpy as np
from icecream import ic
from pathlib import Path
from censo_ext.Tools.utility import is_exist, is_exist_return_bool
# from dataclasses import dataclass


class Anmrrc():
    _instance = None

    def __new__(cls, *args, **kwargs) -> Self:
        '''Singleton Pattern'''
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, DirFileName: Path) -> None:
        lines: list[str] = open(DirFileName, "r").readlines()

        # Dict of Atomic Numbers and Atomic label
        self.element: dict[int, str] = {
            1: 'H', 6: 'C', 9: 'F', 14: 'Si', 15: 'P'}
        # List of Acid Atom (like NH OH) No show in spectra
        self.acid_atoms_NoShow: list = []
        # [atomic number] [calculated shielding valule of the reference molecule] [experimental shift] [active or not]
        self.anmrrc: list[list[float]] = []
        self.Active: list[str] = []
        self.third_line: str = lines[2].rstrip()
        for idx, x in enumerate(lines[0].split()):
            if x != "XH":
                self.acid_atoms_NoShow.append(int(x))
            else:
                break

        second_line: list[str] = lines[1].split()
        self.mf: float = float(second_line[4])             # mf      : nmr frequency           # nopep8
        self.lw: float = float(second_line[6])             # lw      : lines of width          # nopep8
        self.Temp: float = float(second_line[12])          # Temp    : Temperature (K)         # nopep8
        self.JCoups: bool = bool(second_line[8])           # JCoups  : bool of JCoups ONOFF    # nopep8
        self.SParams: bool = bool(second_line[10])         # Sparams : bool of SParams ONOFF   # nopep8

        # 3 lines of .anmrrc Parameters
        # self.third_line: str = lines[2].rstrip()

        for idx, x in enumerate(lines):
            if idx >= 3:
                self.anmrrc.append([int(x.split()[0]), float(
                    x.split()[1]), float(x.split()[2]), int(x.split()[3])])

        # Set the active species based on anmrrc entries
        for idx, x in enumerate(self.anmrrc):
            if x[3] == 1:
                self.Active.append(self.element[int(x[0])])

    def __repr__(self) -> str:
        Result: str = ""
        for x in (self.acid_atoms_NoShow):
            Result += f'{x} '
        Result += f'XH acid atoms\n'
        Result += f'ENSO qm= ORCA mf= {self.mf} lw= {self.lw}  J='
        Result += f" on" if self.JCoups else f" off"
        Result += f" S="
        Result += f" on" if self.SParams else f" off"
        Result += f" T= {self.Temp}\n"
        Result += f"{self.third_line}\n"

        for x in self.anmrrc:
            Result += f"{x[0]:<d}  {x[1]:<9.2f} {x[2]:<6.1f} {x[3]:>2d}\n"
        return Result

    def get_idx1_acid_atoms_NoShow_RemoveH(self, DirFileName: Path) -> list[int]:
        """
            Return the list of be removed hydrogen of acid_atoms_NoShow in .anmrrc setting

        Args:
            fileName (str, optional): _description_. Defaults to "crest_conformers.xyz".
        """

        from censo_ext.Tools.Parameter import ELEMENT_NAMES
        acid_atoms_NoShow: list = [ELEMENT_NAMES[i]
                                   for i in self.acid_atoms_NoShow]
        from censo_ext.Tools.ml4nmr import read_mol_neighbors
        mol, neighbors = read_mol_neighbors(DirFileName)
        acid_atoms_NoShowRemove: list = []
        for idx, x in enumerate(mol):
            for y in acid_atoms_NoShow:
                if x.symbol == y:  # type: ignore
                    acid_atoms_NoShowRemove.append(idx+1)
        NoShow_Remove_Group: np.ndarray = np.array([], dtype=int)
        for idx, x in enumerate(acid_atoms_NoShowRemove):
            NoShow_Remove_Group = np.concatenate(
                (NoShow_Remove_Group, neighbors[x]), axis=None)
        idx_H_atom = [idx+1 for idx,
                      i in enumerate(mol) if i.symbol == "H"]  # type: ignore
        return [i for i in NoShow_Remove_Group if i in idx_H_atom]

    def get_Reference_anmrrc(self) -> float:
        """
            Retrieve the reference shielding value from the .anmrrc file.

            Returns:
                float: The reference shielding value.

            Raises:
                ValueError: If no active species is found in the .anmrrc file.
        """
        reference = None
        for x in self.anmrrc:
            if x[3] == 1:
                reference = x[1]

        if reference is None:
            raise ValueError("No reference in your .anmrrc file")
            exit(0)
        return reference


class Anmr():
    _instance = None

    def __new__(cls, *args, **kwargs) -> Self:
        '''Singleton Pattern'''
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, Directory_str: Path = Path(".")) -> None:
        self.__Directory: Path = Directory_str
        self.enso: np.ndarray                   # anmr_enso
        self.anmrJ: np.ndarray                  # JCoup     of anmr.out generated from anmr
        self.anmrS: list = []                   # Shielding of anmr.out generated from anmr
        self.orcaSJ: list[OrcaSJ] = []          # directory of orcaSJ
        self.avg_orcaSJ = OrcaSJ()          #
        self.nucinfo: list = []                 # anmr_nucinfo

    def get_Directory(self) -> Path:
        return self.__Directory

    def get_Anmr_Active(self) -> list[str]:
        return self.__AnmrParams.Active

    def get_idx1_acid_atoms_NoShow_RemoveH(self, DirFileName=Path("crest_conformers.xyz")) -> list[int]:
        return self.__AnmrParams.get_idx1_acid_atoms_NoShow_RemoveH(DirFileName)

    def method_read_anmrrc(self, fileName=Path(".anmrrc")) -> None:
        '''Read .anmrrc setting file from censo '''
        DirFileName: Path = self.__Directory / fileName
        is_exist(DirFileName)
        self.__AnmrParams: Anmrrc = Anmrrc(DirFileName)
        # self.__AnmrParams: Anmrrc = Anmrrc(lines)

    def method_print_anmrrc(self) -> None:
        print(self.__AnmrParams, end="")

    def method_avg_orcaSJ(self) -> None:
        '''Average of orcaJ.out and orcaS.out in /NMR/CONFXX folder'''
        print(" ===== Average of all folder orcaS.out and orcaJ.out =====")

        if len(self.nucinfo) == 0 or len(self.enso) == 0:
            print(" Need to read the anmr_nucinfo and anmr_enso enso")
            print(" Exit the program !!! ")
            ic()
            exit(0)
        else:
            # for Normal of weight of anmr_enso
            weight = self.enso['BW']
            switch = self.enso['ONOFF']
            idx_CONF = self.enso['CONF'].astype(int)

            if (np.sum(switch) == 0):
                print("anmr_enso: Table - ONOFF is Zero ")
                ic()
                exit(0)

            weight = weight*switch
            weight = weight / sum(weight)
            normal_weight = dict(zip(idx_CONF, weight))

            import copy
            Active_orcaSJ: list = []
            for idx, x in enumerate(self.orcaSJ):
                if x.CONFSerialNums in idx_CONF:
                    Active_orcaSJ.append(idx)

            # orcaSParams and orcaJCoups using weighting to calculate and
            # save to Average_orcaSJ
            self.avg_orcaSJ: OrcaSJ = OrcaSJ()
            self.avg_orcaSJ.idxAtoms = self.orcaSJ[0].idxAtoms

            for idx in self.orcaSJ[0].SParams.keys():
                self.avg_orcaSJ.SParams[idx] = 0.0

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                y_idx = list(x.SParams.keys())
                y = list(x.SParams.values())
                for idz, z in zip(y_idx, np.array(y) * normal_weight[x.CONFSerialNums]):
                    self.avg_orcaSJ.SParams[idz] += z

            for key, value in self.avg_orcaSJ.SParams.items():
                self.avg_orcaSJ.SParams[key] = value - \
                    self.__AnmrParams.get_Reference_anmrrc()

            nNumbers: int = np.shape(self.orcaSJ[0].JCoups[0])[0]
            self.avg_orcaSJ.JCoups = np.zeros((nNumbers, nNumbers))

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                y = x.JCoups
                self.avg_orcaSJ.JCoups += y * \
                    normal_weight[x.CONFSerialNums]

            ic(self.avg_orcaSJ.SParams)
            ic(self.avg_orcaSJ.JCoups)

        print(" ===== Finished the Average of all folder orcaS.out and orcaJ.out =====")

    def method_update_equiv_orcaSJ(self) -> None:
        '''
        update the equivalent of SParams and JCoups in /NMR/CONFXX folder
        '''
        print(" Replace the equivalent of Sparams and JCoups")

        if len(self.nucinfo) == 0 or len(self.enso) == 0:
            print(" Need to read the anmr_nucinfo and anmr_enso enso")
            print(" Exit the program !!! ")
            ic()
            exit(1)
        else:
            print("===== Update the equivalent of SParams and JCoups =====")

            List_tmp: list[int] = []
            for idx, x in self.orcaSJ[0].idxAtoms.items():
                List_tmp.append(idx)
            AtomsKeep: np.ndarray = np.array(List_tmp)
            AtomsEqvKeep: np.ndarray = AtomsKeep.copy()

            ic(AtomsEqvKeep)
            ic(self.nucinfo)

            for idx, x in enumerate(AtomsEqvKeep):
                if (self.nucinfo[1][x-1][1]) != 1:
                    if x != min(self.nucinfo[1][x-1][2]):
                        AtomsEqvKeep[idx] = 0
            AtomsEqvKeep = AtomsEqvKeep[np.nonzero(AtomsEqvKeep)]

            # Calculation the average ppm of Equivalent Atom and Replace the old ppm
            for orcaSJ in self.orcaSJ:
                for idx, x in enumerate(AtomsEqvKeep):
                    if (self.nucinfo[0][x-1][1] != 1):
                        ppm: list[float] = []
                        for idy, y in enumerate(self.nucinfo[0][x-1][2]):
                            ppm.append(orcaSJ.SParams[y])
                        average: float = sum(ppm)/len(ppm)
                        for y in self.nucinfo[0][x-1][2]:
                            orcaSJ.SParams[y] = average

            # for Equivalent Atom  of orcaJCoups
            # Calculation the average JCoups of Equivalent JCoups and Replace the old JCoups
            for orcaSJ in self.orcaSJ:
                for idx, x in enumerate(AtomsEqvKeep):
                    if (self.nucinfo[1][x-1][1] != 1):
                        for idy, y in enumerate(AtomsKeep):
                            JCoups: list[float] = []
                            average: float = 0
                            for z in (self.nucinfo[1][x-1][2]):
                                JCoups.append(orcaSJ.JCoups[list(
                                    AtomsKeep).index(y)][list(AtomsKeep).index(z)])
                                average: float = sum(JCoups)/len(JCoups)

                            for z in (self.nucinfo[1][x-1][2]):
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    y)][list(AtomsKeep).index(z)] = average
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    z)][list(AtomsKeep).index(y)] = average

                    if (self.nucinfo[1][x-1][1] > 2):
                        for k in (self.nucinfo[1][x-1][2]):
                            for l in (self.nucinfo[1][x-1][2]):
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    k)][list(AtomsKeep).index(l)] = 0
                                orcaSJ.JCoups[list(AtomsKeep).index(
                                    l)][list(AtomsKeep).index(k)] = 0
                for idx, x in enumerate(orcaSJ.JCoups):
                    orcaSJ.JCoups[idx][idx] = 0

            AtomsDelete: list[int] = list(
                set(AtomsKeep).difference(set(AtomsEqvKeep)))
            AtomsDelete.sort()

            # Delete Equivalent Atoms of idxAtoms
            for orcaSJ in self.orcaSJ:
                for x in AtomsDelete:
                    del orcaSJ.idxAtoms[x]

            # Delete Equivalent Atoms of orcaSParams
            for orcaSJ in self.orcaSJ:
                for x in AtomsDelete:
                    del orcaSJ.SParams[x]

            # Delete Equivalent Atom of orcaJCoups
            AtomsDelete2idx0: dict[int, int] = {}
            for idx, x in enumerate(AtomsKeep):
                if x in AtomsDelete:
                    AtomsDelete2idx0[x] = idx
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
                        if not idy0 in idx0_AtomsDelete:
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

        for idx in range(len(dirNames)):
            file_orcaS: Path = Dir / Path(dirNames[idx] + "/NMR/orcaS.out")  # nopep8
            file_orcaJ: Path = Dir / Path(dirNames[idx] + "/NMR/orcaJ.out")  # nopep8
            if (os.path.exists(file_orcaS) == True and os.path.exists(file_orcaJ) == True):
                print(str(idx)+"  :  "+str(file_orcaS))
                print(str(idx)+"  :  "+str(file_orcaJ))

                iter: OrcaSJ = OrcaSJ()
                iter.CONFSerialNums = int(dirNames[idx].replace('CONF', ''))
                if (iter.method_read_orcaS(file=file_orcaS) == False):
                    print("Something wrong in your orcaS.out")
                if (iter.method_read_orcaJ(file=file_orcaJ) == False):
                    print("Something wrong in your orcaJ.out")
                self.orcaSJ.append(iter)

    def get_avg_orcaSJ_Exist(self) -> bool:

        fileName_Av_orcaS: Path = self.__Directory / \
            Path("Average/NMR/orcaS.out")
        fileName_Av_orcaJ: Path = self.__Directory / \
            Path("Average/NMR/orcaJ.out")
        fileName_Av_orcaAtoms: Path = self.__Directory / \
            Path("Average/NMR/orcaA.out")

        if is_exist_return_bool(fileName_Av_orcaS) and is_exist_return_bool(fileName_Av_orcaAtoms) and is_exist_return_bool(fileName_Av_orcaJ):
            return True
        else:
            return False

    def method_load_avg_orcaSJ(self, args) -> bool:

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
                self.avg_orcaSJ.idxAtoms = json.loads(
                    f.read(), object_pairs_hook=jsonKeys2int)

            from censo_ext.Tools.utility import load_dict_orcaS
            self.avg_orcaSJ.SParams = load_dict_orcaS(
                file_avg_orcaS)
            self.avg_orcaSJ.JCoups = np.loadtxt(file_avg_orcaJ)
            return True
        else:
            return False

    def method_save_adjust_avg_orcaS(self) -> None:
        raise NotImplementedError
        # fileName_Average_orcaS = str(
        #    self.__Directory / Path("Average/NMR/orcaS.out"))
        # from censo_ext.Tools.utility import Save_Dict_orcaS
        # Save_Dict_orcaS(fileName_Average_orcaS,
        #                self.Average_orcaSJ.orcaSParams)

    def method_save_avg_orcaSJ(self) -> None:

        avg_orcaS: Path = self.__Directory / Path("Average/NMR/orcaS.out")      # nopep8
        avg_orcaJ: Path = self.__Directory / Path("Average/NMR/orcaJ.out")      # nopep8
        avg_orcaAtoms: Path = self.__Directory / Path("Average/NMR/orcaA.out")  # nopep8

        Path(self.__Directory / Path("Average/NMR")
             ).mkdir(parents=True, exist_ok=True)

        from censo_ext.Tools.utility import save_dict_orcaS
        save_dict_orcaS(avg_orcaS, self.avg_orcaSJ.SParams)
        import json
        with open(avg_orcaAtoms, 'w') as f:
            f.write(json.dumps(self.avg_orcaSJ.idxAtoms))
        np.savetxt(avg_orcaJ, self.avg_orcaSJ.JCoups, fmt="%10.5f")

    def method_save_folder_orcaSJ(self) -> None:
        raise NotImplementedError("Under Construct")

    def method_read_anmrSJ(self, filename: Path = Path("anmr.out")) -> None:
        '''
        Read the file anmr.out from anmr program
        '''
        filename = self.__Directory / Path(filename)
        is_exist(filename)

        start_idx, end_idx = 0, 0
        DataJ: list[str] = []
        lines: list[str] = open(filename, "r").readlines()

        import re
        bool_first_line = False
        start_shielding_idx: int = 0
        for idx, line in enumerate(lines):
            if re.search(r"MATRIX PRINTED:", line):
                start_idx = idx+1
            if re.search(r"\+\/\-", line) and bool_first_line == False:
                start_shielding_idx = idx
                bool_first_line = True
        if start_shielding_idx != 0:
            nNuclei: int = start_idx - start_shielding_idx - 2
        else:
            print("Something wrong")
            ic()
            exit(1)

        nLines = 0
        for idx in range(int(nNuclei/6)):
            nLines += nNuclei - idx * 6 + 3
        end_idx = start_idx + nLines + nNuclei % 6 + 3 - 1

        for idx, line in enumerate(lines):
            if idx >= start_idx and idx <= end_idx:
                DataJ.append(line.rstrip())

        for idx, line in enumerate(lines):
            if re.search(r"\+\/\-", line):
                tmp: list[float] = [int(i) for i in line.split()[0:3]]
                tmp.append(float(line.split()[3]))
                self.anmrS.append(tmp)

        for idx, line in enumerate(lines):
            if re.search(r"1H resonance frequency", line):
                self.frq = float(line.split()[6])

        ListDataJ: list[str] = DataJ
        nLinesDataJ: int = len(self.anmrS)
        ListDataJDeleteBlank: list = []

        k: int = nLinesDataJ
        i: int = len(ListDataJ)

        # Delete the top serial numbers of JCoup
        while i >= 1:
            del ListDataJ[0:3]
            ListDataJDeleteBlank += ListDataJ[0:k]
            ListDataJ = ListDataJ[k:]
            i, k = i-6, k-6

        # Delete the right serial numbers of JCoup
        ListDataDelete_nAtoms: list[str] = []
        for line in ListDataJDeleteBlank:
            if (len(ListDataDelete_nAtoms) % nLinesDataJ == 0):
                for j in range(6*((int)(len(ListDataDelete_nAtoms)/nLinesDataJ))):
                    ListDataDelete_nAtoms.append("\n")
            ListDataDelete_nAtoms.append(line[6:])

        DataJ_triangle = [""]*(nLinesDataJ)

        # Restruct to one all JCoup table
        for idx, line in enumerate(ListDataDelete_nAtoms):
            DataJ_triangle[idx % nLinesDataJ] = DataJ_triangle[idx %
                                                               nLinesDataJ].rstrip("\n") + " " + line.rstrip("\n")
        self.anmrJ = np.zeros((nLinesDataJ, nLinesDataJ))

        # copy half to other half data on JCoup
        for idx in range(nLinesDataJ):
            for idy in range(idx):
                self.anmrJ[idy][idx] = self.anmrJ[idx][idy] = float(
                    DataJ_triangle[idx].split()[idy])

    def method_print_anmrS(self) -> None:
        print("    #  in coord file  # nucs   delta(ppm)")
        for idx, x in enumerate(self.anmrS):
            print(f"{x[0]:>5d} {x[1]:>9d} {x[2]:>9d} {x[3]:>13.3f}")

    def method_print_anmrJ(self) -> None:
        for idx in range(self.anmrJ[0].size):
            for idy in range(self.anmrJ[0].size):
                print(f'{self.anmrJ[idx][idy]:>10.5f}', end="")
            print("")

    def method_print_nucinfo(self) -> None:
        nAtoms: int = len(self.nucinfo[0])
        print(f"{nAtoms:>12d}")

        for lines in self.nucinfo:
            for x in lines:
                print(f"   {x[0]:d}   {x[1]:d}")
                for y in x[2]:
                    print(f" {y:d}", end="")
                print("")

    def method_read_nucinfo(self, file: Path = Path("anmr_nucinfo")) -> None:
        file = self.__Directory / Path(file)
        is_exist(file)

        lines: list[str] = open(file, "r").readlines()
        nAtoms: int = int(lines[0].split()[0])
        del lines[0]
        page: list = []
        for idx, x in enumerate(lines):
            x = x.rstrip()
            if (idx % 2) == 0:
                tmp = []
                tmp.append(int(x.split()[0]))
                tmp.append(int(x.split()[1]))
            else:
                int_tmp = []
                for idy, y in enumerate(x.split()):
                    int_tmp.append(int(y))
                tmp.append(int_tmp)         # type: ignore
                page.append(tmp)            # type: ignore
            if ((idx+1) % (nAtoms*2)) == 0:
                self.nucinfo.append(page)
                page = []

    def method_create_enso(self, in_np: np.ndarray) -> None:
        self.enso = in_np
        if len(self.enso[0]) != 8:
            print("something wrong in your anmr_enso file")
            ic()
            exit(1)

    def method_read_enso(self, file: Path = Path("anmr_enso")) -> None:
        ''' anmr_enso :  8 columns '''
        ''' ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi  '''

        file = self.__Directory / Path(file)
        is_exist(file)

        self.enso = np.genfromtxt(file, names=True)

        if len(self.enso[0]) != 8:
            print("something wrong in your anmr_enso file")
            ic()
            exit(1)

    def method_print_enso(self) -> None:
        print("ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi     ")
        for Enso in self.enso:
            print(f'{int(Enso[0]):<1d}      {int(Enso[1]):<4d} {int(Enso[2]):<4d} {Enso[3]:>6.4f} {Enso[4]: > 11.7f} {Enso[5]: > 10.7f} {Enso[6]: > 10.7f} {Enso[7]:>4.3f}')  # nopep8

    def method_save_enso(self, file: Path = Path("anmr_enso.new")) -> None:
        DirFileName: Path = self.__Directory / Path(file)
        original_stdout = sys.stdout
        with open(DirFileName, "w") as f:
            sys.stdout = f
            self.method_print_enso()
        sys.stdout = original_stdout


class OrcaSJ():
    def __init__(self) -> None:
        self.JCoups: np.ndarray = np.array([])
        self.SParams: dict[int, float] = {}
        self.Anisotropy: dict[int, float] = {}
        self.CONFSerialNums: int = 0
        self.idxAtoms: dict[int, str] = {}

    def method_read_orcaJ(self, file: Path = Path("orcaJ.out")) -> bool:
        ''' Read the file orcaJ.out in censo program '''
        print(" method_read_orcaJ", file)
        is_exist(file)

        start_idx, end_idx = 0, 0
        DataJ: list = []
        lines: list[str] = open(file, "r").readlines()
        import re
        nVersion: str = ""

        for idx, line in enumerate(lines):
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
                    # end_idx = idx-2

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
            print("    exit and close the program !!! ")
            ic()
            exit(0)

        for idx, line in enumerate(lines):
            if idx >= start_idx and idx <= end_idx:
                DataJ.append(line.rstrip())

        for idx, x in enumerate(DataJ):
            DataJ[idx] = x.split()[2:]
        nums = 1
        while (nums >= 1):
            if ((int((nums-1)/6)+1)*(nums+1) == len(DataJ)):
                break
            nums = nums+1
        nAtomDataJ = nums

        for idx, x in enumerate(DataJ):
            if (idx > nAtomDataJ):
                DataJ[idx % (nAtomDataJ+1)] = DataJ[idx % (nAtomDataJ+1)]+x

        del DataJ[0]
        del DataJ[nAtomDataJ:]
        self.JCoups = np.array(DataJ, dtype=float)
        return True

    def method_read_orcaS(self, file: Path = Path("orcaS.out")) -> bool:
        ''' Read the file orcaS.out in censo program '''
        print(" method_read_orcaS", file)
        is_exist(file)

        start_idx, end_idx = 0, 0
        DataS: list = []
        lines: list[str] = open(file, "r").readlines()
        import re
        nVersion: str = ""
        nNuclei = 0

        for idx, line in enumerate(lines):
            if re.search(r"Program Version", line):
                nVersion: str = line

        if int(nVersion.split()[2][0]) == 5:
            for idx, line in enumerate(lines):
                if re.search(r"Number of nuclei for epr/nmr", line):
                    nNuclei = int(line.split()[-1])
                if re.search(r"CHEMICAL SHIELDING SUMMARY", line):
                    start_idx = idx + 6
            end_idx = start_idx + nNuclei - 1
            if end_idx == 0 or start_idx == 0 or nNuclei == 0:
                return False
        elif int(nVersion.split()[2][0]) == 6:
            for idx, line in enumerate(lines):
                if re.search(r"Maximum memory used throughout the entire PROP", line):
                    end_idx = idx - 5
                if re.search(r"CHEMICAL SHIELDING SUMMARY", line):
                    start_idx = idx + 6
        else:
            print("This program is not work with before orca 5.0 ")
            ic()
            exit(0)

        for idx, line in enumerate(lines):
            if idx >= start_idx and idx <= end_idx:
                DataS.append(line.rstrip())
        self.idxAtoms, self.Anisotropy, self.SParams = {}, {}, {}
        for idx, x in enumerate(DataS):
            self.idxAtoms[int(x.split()[0])+1] = str(x.split()[1])
            self.SParams[int(x.split()[0])+1] = float(x.split()[2])
            self.Anisotropy[int(x.split()[0])+1] = float(x.split()[3])
        return True

    def method_save_orcaS(self) -> list:
        raise NotImplementedError
        # if len(self.idxAtoms) == len(self.SParams):
        #    lines: list = []
        #    lines.append("  Nucleus  Element    Isotropic     Anisotropy\n")
        #    lines.append("  -------  -------  ------------   ------------\n")
        #
        #    for key, value in self.idxAtoms.items():
        #        str1: str = f"  \
        #            {key-1: > 5d}{value: > 8s}{self.SParams[key]: > 15.3f}{self.Anisotropy[key]: > 15.3f}\n"
        #        lines.append(str1)
        #    return lines
        #
        # else:
        #    print("your orcaJ and orcaS is not fit each other")
        #    print("    exit and close the program !!! ")
        #    ic()
        #    exit(0)

    def method_save_orcaJ(self) -> list:
        raise NotImplementedError
        if len(self.idxAtoms) == len(self.JCoups[0]):
            lines: list = []
            list_idxAtoms: list = list(self.idxAtoms.keys())

            for idx in range(len(self.idxAtoms)):
                for idy in range(idx+1, len(self.idxAtoms), 1):
                    # print(idx, idy)
                    str1: str = f" NUCLEUS A = {(self.idxAtoms[list_idxAtoms[idx]]): > s}\
                        {(list_idxAtoms[idx]-1): > 5d} NUCLEUS B = {self.idxAtoms[list_idxAtoms[idy]]: > s}{list_idxAtoms[idy]-1: > 5d}"
                    # " NUCLEUS A = H    4 NUCLEUS B = H    5"
                    str2: str = f" Total            0.000            0.000            0.000  iso = \
                        {self.JCoups[idx][idy]: > 13.5f}"
                    # " Total            4.506            4.506            4.506  iso=       4.506"
                    lines.append(str1)
                    lines.append(str2)
            return lines
        else:
            print("your orcaJ and orcaS is not fit each other")
            print("    exit and close the program !!! ")
            ic()
            exit(0)

    def method_print_orcaS(self) -> None:
        if len(self.idxAtoms) == len(self.SParams):
            print("Nucleus  Element   Anisotropy")
            for idx, Atom in self.idxAtoms.items():
                print(f'{idx:>5d}', f'{Atom:>8s}', end="")
                print(f'{self.SParams[idx]:>15.3f}')
        else:
            print("your orcaJ and orcaS is not fit each other")
            print("    exit and close the program !!! ")
            ic()
            exit(0)

    def method_print_orcaJ(self) -> None:
        for idx in range(self.JCoups[0].size):
            for idy in range(self.JCoups[0].size):
                print(f'{(self.JCoups[idx][idy]):>8.3f}', end="")
            print("")


class CensoDat():

    def __init__(self, file: Path = Path("anmr.dat")) -> None:
        is_exist(file)
        self.__fileName: Path = Path(file)
        self.__dat: np.ndarray = np.genfromtxt(file)

    def __len__(self) -> int:
        return len(self.__dat)

    def __sub__(self, other: Self) -> CensoDat:
        if np.array_equal(self.__dat[:, 0], other.__dat[:, 0]):
            self.__dat[:, 1] -= other.__dat[:, 1]
        return self

    def __repr__(self) -> str:
        Result: str = ""
        for x in self.__dat:
            Result = Result + f'{x[0]:>12.6f}  {x[1]:>12.6e}\n'
        return Result

    def method_save_dat(self) -> None:
        original_stdout = sys.stdout
        with open(self.__fileName, "w") as f:
            sys.stdout = f
            print(self, end="")
        sys.stdout = original_stdout

    def method_normalize_dat(self, start: float = -5.0, end: float = 15.0, dpi: int = 10000, highest: int = 10000) -> None:

        ppm_least: float = self.__dat[:, 0][-1]
        if ppm_least > 50 and end < 50:
            start, end, dpi = -20, 240, 500

        if len(self) != 0:

            res: np.ndarray = self.__dat
            res = np.insert(res, 0, [start, 0.0], axis=0)
            res = np.vstack((res, [end, 0.0]))

            from scipy import interpolate
            f = interpolate.interp1d(res[:, 0], res[:, 1])

            xnew: np.ndarray = np.linspace(start, end, int(end-start)*dpi+1)
            ynew: np.ndarray = f(xnew)

            res_new: np.ndarray = np.vstack((xnew, ynew))
            res_new[1] = res_new[1] / np.max(res_new[1]) * highest
            self.__dat = res_new.T

    def method_subtract_dat(self, other: Self) -> CensoDat:
        if np.array_equal(self.__dat[:, 0], other.__dat[:, 0]):
            self.__dat[:, 1] -= other.__dat[:, 1]
        else:
            print("two dat file is not the same scale")
            ic()
            exit(0)
        return self

    def set_fileName(self, file: Path) -> None:
        self.__fileName = Path(file)

    def get_fileName(self) -> Path:
        return self.__fileName

    def get_Dat(self) -> np.ndarray:
        return self.__dat
