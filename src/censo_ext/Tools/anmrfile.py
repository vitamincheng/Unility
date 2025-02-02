#!/usr/bin/env python3
#  Module anmr/ censo           [08.27.2024] vitamin.cheng@gmail.com
from __future__ import annotations
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
        # self.__AnmrParams: Anmrrc = Anmrrc(lines)

        # Dict of Atomic Numbers and Atomic label
        self.element: dict[int, str] = {
            1: 'H', 6: 'C', 9: 'F', 14: 'Si', 15: 'P'}
        # List of Acid Atom (like NH OH) No show in spectra
        self.AcidAtomsNoShow: list = []
        # [atomic number] [calculated shielding valule of the reference molecule] [experimental shift] [active or not]
        self.anmrrc: list[list[float]] = []
        self.Active: list[str] = []
        self.Thirdline: str = lines[2].rstrip()
        for idx, x in enumerate(lines[0].split()):
            if x != "XH":
                self.AcidAtomsNoShow.append(int(x))
            else:
                break

        second_line: list[str] = lines[1].split()
        self.mf: float = float(second_line[4])             # mf      : nmr frequency           # nopep8
        self.lw: float = float(second_line[6])             # lw      : lines of width          # nopep8
        self.Temp: float = float(second_line[12])          # Temp    : Temperature (K)         # nopep8
        self.JCoups: bool = bool(second_line[8])           # JCoups  : bool of JCoups ONOFF    # nopep8
        self.SParams: bool = bool(second_line[10])         # Sparams : bool of SParams ONOFF   # nopep8

        # 3 lines of .anmrrc Parameters
        # self.Thirdline: str = lines[2].rstrip()

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
        for x in (self.AcidAtomsNoShow):
            Result += f'{x} '
        Result += f'XH acid atoms\n'
        Result += f'ENSO qm= ORCA mf= {self.mf} lw= {self.lw}  J='
        Result += f" on" if self.JCoups else f" off"
        Result += f" S="
        Result += f" on" if self.SParams else f" off"
        Result += f" T= {self.Temp}\n"
        Result += f"{self.Thirdline}\n"

        for x in self.anmrrc:
            Result += f"{x[0]:<d}  {x[1]:<9.2f} {x[2]:<6.1f} {x[3]:>2d}\n"
        return Result

    def get_list_idx1AcidAtomsNoShowRemoveH(self, DirFileName: Path) -> list[int]:
        """
            Return the list of be removed hydrogen of AcidAtomsNoShow in .anmrrc setting

        Args:
            fileName (str, optional): _description_. Defaults to "crest_conformers.xyz".
        """

        from censo_ext.Tools.Parameter import ELEMENT_NAMES
        AcidAtomsNoShow: list = [ELEMENT_NAMES[i]
                                 for i in self.AcidAtomsNoShow]
        from censo_ext.Tools.ml4nmr import read_mol_neighbors
        mol, neighbors = read_mol_neighbors(DirFileName)
        AcidAtomsNoShowRemove: list = []
        for idx, x in enumerate(mol):
            for y in AcidAtomsNoShow:
                if x.symbol == y:  # type: ignore
                    AcidAtomsNoShowRemove.append(idx+1)
        AcidAtomsNoShowRemoveGroup: np.ndarray = np.array([], dtype=int)
        for idx, x in enumerate(AcidAtomsNoShowRemove):
            AcidAtomsNoShowRemoveGroup = np.concatenate(
                (AcidAtomsNoShowRemoveGroup, neighbors[x]), axis=None)
        idx_H_atom = [idx+1 for idx,
                      i in enumerate(mol) if i.symbol == "H"]  # type: ignore
        return [i for i in AcidAtomsNoShowRemoveGroup if i in idx_H_atom]

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
        self.Average_orcaSJ = OrcaSJ()          #
        self.nucinfo: list = []                 # anmr_nucinfo

    def get_Directory(self) -> Path:
        return self.__Directory

    def get_Anmr_Active(self) -> list[str]:
        return self.__AnmrParams.Active

    def get_list_idx1AcidAtomsNoShowRemoveH(self, DirFileName=Path("crest_conformers.xyz")) -> list[int]:
        return self.__AnmrParams.get_list_idx1AcidAtomsNoShowRemoveH(DirFileName)

    def method_read_anmrrc(self, fileName=Path(".anmrrc")) -> None:
        '''Read .anmrrc setting file from censo '''
        DirFileName: Path = self.__Directory / fileName
        is_exist(DirFileName)
        self.__AnmrParams: Anmrrc = Anmrrc(DirFileName)
        # self.__AnmrParams: Anmrrc = Anmrrc(lines)

    def method_print_anmrrc(self) -> None:
        print(self.__AnmrParams, end="")

    def method_average_orcaSJ(self) -> None:
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

            # ic(switch)
            weight = weight*switch
            weight = weight / sum(weight)
            normal_weight = dict(zip(idx_CONF, weight))
            # ic(normal_weight)

            import copy
            Active_orcaSJ: list = []
            for idx, x in enumerate(self.orcaSJ):
                if x.CONFSerialNums in idx_CONF:
                    Active_orcaSJ.append(idx)
            # ic(Active_orcaSJ)

            # orcaSParams and orcaJCoups using weighting to calculate and
            # save to Average_orcaSJ
            self.Average_orcaSJ: OrcaSJ = OrcaSJ()
            self.Average_orcaSJ.idxAtoms = self.orcaSJ[0].idxAtoms

            for idx in self.orcaSJ[0].orcaSParams.keys():
                self.Average_orcaSJ.orcaSParams[idx] = 0.0

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                y_idx = list(x.orcaSParams.keys())
                y = list(x.orcaSParams.values())
                for idz, z in zip(y_idx, np.array(y) * normal_weight[x.CONFSerialNums]):
                    self.Average_orcaSJ.orcaSParams[idz] += z

            for key, value in self.Average_orcaSJ.orcaSParams.items():
                self.Average_orcaSJ.orcaSParams[key] = value - \
                    self.__AnmrParams.get_Reference_anmrrc()

            nNumbers: int = np.shape(self.orcaSJ[0].orcaJCoups[0])[0]
            self.Average_orcaSJ.orcaJCoups = np.zeros((nNumbers, nNumbers))

            for x in np.array(self.orcaSJ)[Active_orcaSJ]:
                y = x.orcaJCoups
                self.Average_orcaSJ.orcaJCoups += y * \
                    normal_weight[x.CONFSerialNums]

            ic(self.Average_orcaSJ.orcaSParams)
            ic(self.Average_orcaSJ.orcaJCoups)

        print(" ===== Finished the Average of all folder orcaS.out and orcaJ.out =====")

    def method_update_equivalent_orcaSJ(self) -> None:
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
            ic(AtomsKeep)
            ic(AtomsEqvKeep)

            # Calculation the average ppm of Equivalent Atom and Replace the old ppm
            for orcaSJ in self.orcaSJ:
                for idx, x in enumerate(AtomsEqvKeep):
                    if (self.nucinfo[0][x-1][1] != 1):
                        ppm: list[float] = []
                        for idy, y in enumerate(self.nucinfo[0][x-1][2]):
                            ppm.append(orcaSJ.orcaSParams[y])
                        average: float = sum(ppm)/len(ppm)
                        for y in self.nucinfo[0][x-1][2]:
                            orcaSJ.orcaSParams[y] = average

            # ic(self.orcaSJ[2].orcaSParams)
            # ic(self.nucinfo[0])

            # for Equivalent Atom  of orcaJCoups
            # Calculation the average JCoups of Equivalent JCoups and Replace the old JCoups
            # print(" Average of orcaJCoups : \n", self.Average_orcaSJ.orcaJCoups)

            for orcaSJ in self.orcaSJ:
                # print(self.orcaSJ[num].orcaJCoups)
                for idx, x in enumerate(AtomsEqvKeep):
                    if (self.nucinfo[1][x-1][1] != 1):
                        for idy, y in enumerate(AtomsKeep):
                            # print(x, y, list(AtomsKeep).index(y))
                            JCoups: list[float] = []
                            average: float = 0
                            for z in (self.nucinfo[1][x-1][2]):
                                JCoups.append(orcaSJ.orcaJCoups[list(
                                    AtomsKeep).index(y)][list(AtomsKeep).index(z)])
                                average: float = sum(JCoups)/len(JCoups)

                            for z in (self.nucinfo[1][x-1][2]):
                                orcaSJ.orcaJCoups[list(AtomsKeep).index(
                                    y)][list(AtomsKeep).index(z)] = average
                                orcaSJ.orcaJCoups[list(AtomsKeep).index(
                                    z)][list(AtomsKeep).index(y)] = average

                    if (self.nucinfo[1][x-1][1] > 2):
                        for k in (self.nucinfo[1][x-1][2]):
                            for l in (self.nucinfo[1][x-1][2]):
                                orcaSJ.orcaJCoups[list(AtomsKeep).index(
                                    k)][list(AtomsKeep).index(l)] = 0
                                orcaSJ.orcaJCoups[list(AtomsKeep).index(
                                    l)][list(AtomsKeep).index(k)] = 0

                for idx, x in enumerate(orcaSJ.orcaJCoups):
                    orcaSJ.orcaJCoups[idx][idx] = 0
                # print(self.orcaSJ[num].orcaJCoups)
            # for orcaSJ in self.orcaSJ:
            #    ic(orcaSJ.CONFSerialNums)
            #    ic(orcaSJ.orcaSParams)
            #    ic(orcaSJ.orcaJCoups)

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
                    del orcaSJ.orcaSParams[x]

            # Delete Equivalent Atom of orcaJCoups
            ic(AtomsDelete)
            ic(AtomsKeep)
            AtomsDelete2idx0: dict[int, int] = {}
            # AtomsDelete2 = AtomsDelete
            # AtomsDelete2.reverse()
            for idx, x in enumerate(AtomsKeep):
                if x in AtomsDelete:
                    AtomsDelete2idx0[x] = idx
            ic(AtomsDelete2idx0)
            list_AtomsDelete: list[int] = [
                x for x in AtomsDelete2idx0.values()]
            list_AtomsDelete.reverse()
            ic(list_AtomsDelete)

            for orcaSJ in self.orcaSJ:
                for x in list_AtomsDelete:
                    orcaSJ.orcaJCoups = np.delete(orcaSJ.orcaJCoups, x, 0)
                    orcaSJ.orcaJCoups = np.delete(orcaSJ.orcaJCoups, x, 1)

            List_idx1AcidAtomsNoShowRemoveH: list[int] = self.__AnmrParams.get_list_idx1AcidAtomsNoShowRemoveH(
                self.__Directory / Path("crest_conformers.xyz"))

            # Delete orcaSJ SParams in AcidAtomsNoShow
            idx0_AtomsDelete: list[int] = []
            for orcaSJ in self.orcaSJ:
                for idy0, SParam in enumerate(orcaSJ.orcaSParams.copy()):
                    if SParam in List_idx1AcidAtomsNoShowRemoveH:
                        if not idy0 in idx0_AtomsDelete:
                            idx0_AtomsDelete.append(idy0)
                        del orcaSJ.orcaSParams[SParam]

            # Delete orcaSJ orcaJCoups in AcidAtomsNoShow
            ic(idx0_AtomsDelete)
            idx0_AtomsDelete.sort
            idx0_AtomsDelete.reverse()

            for orcaSJ in self.orcaSJ:
                for x in idx0_AtomsDelete:
                    orcaSJ.orcaJCoups = np.delete(orcaSJ.orcaJCoups, x, 0)
                    orcaSJ.orcaJCoups = np.delete(orcaSJ.orcaJCoups, x, 1)

            # for orcaSJ in self.orcaSJ:
            #    ic(orcaSJ.CONFSerialNums)
            #    ic(orcaSJ.orcaSParams)
            #    ic(orcaSJ.orcaJCoups)

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
            filename_orcaS: Path = Dir / Path(dirNames[idx] + "/NMR/orcaS.out")  # nopep8
            filename_orcaJ: Path = Dir / Path(dirNames[idx] + "/NMR/orcaJ.out")  # nopep8
            if (os.path.exists(filename_orcaS) == True and os.path.exists(filename_orcaJ) == True):
                print(str(idx)+"  :  "+str(filename_orcaS))
                print(str(idx)+"  :  "+str(filename_orcaJ))

                iter: OrcaSJ = OrcaSJ()
                iter.CONFSerialNums = int(dirNames[idx].replace('CONF', ''))
                if (iter.method_read_orcaS(filename=filename_orcaS) == False):
                    print("Something wrong in your orcaS.out")
                if (iter.method_read_orcaJ(filename=filename_orcaJ) == False):
                    print("Something wrong in your orcaJ.out")
                self.orcaSJ.append(iter)

                # ic(iter.orcaSParams)
                # ic(iter.orcaJCoups)

    def get_average_orcaSJ_Exist(self) -> bool:

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

    def method_load_average_orcaSJ(self, args) -> bool:

        from censo_ext.Tools.utility import jsonKeys2int
        if self.get_average_orcaSJ_Exist():
            if args.bobyqa:
                filename_Av_orcaS: Path = self.__Directory / \
                    Path("Average/NMR/orcaS-BOBYQA.out")
            else:
                filename_Av_orcaS: Path = self.__Directory / \
                    Path("Average/NMR/orcaS.out")
            filename_Av_orcaJ: Path = self.__Directory / \
                Path("Average/NMR/orcaJ.out")
            filename_Av_orcaAtoms: Path = self.__Directory / \
                Path("Average/NMR/orcaA.out")

            import json
            with open(filename_Av_orcaAtoms) as f:
                self.Average_orcaSJ.idxAtoms = json.loads(
                    f.read(), object_pairs_hook=jsonKeys2int)

            from censo_ext.Tools.utility import load_dict_orcaS
            self.Average_orcaSJ.orcaSParams = load_dict_orcaS(
                filename_Av_orcaS)
            self.Average_orcaSJ.orcaJCoups = np.loadtxt(filename_Av_orcaJ)
            return True
        else:
            return False

    def method_save_adjust_average_orcaS(self) -> None:
        raise NotImplementedError
        # fileName_Average_orcaS = str(
        #    self.__Directory / Path("Average/NMR/orcaS.out"))
        # from censo_ext.Tools.utility import Save_Dict_orcaS
        # Save_Dict_orcaS(fileName_Average_orcaS,
        #                self.Average_orcaSJ.orcaSParams)

    def method_save_average_orcaSJ(self) -> None:

        Av_orcaS: Path = self.__Directory / Path("Average/NMR/orcaS.out")      # nopep8
        Av_orcaJ: Path = self.__Directory / Path("Average/NMR/orcaJ.out")      # nopep8
        Av_orcaAtoms: Path = self.__Directory / Path("Average/NMR/orcaA.out")  # nopep8

        Path(self.__Directory / Path("Average/NMR")
             ).mkdir(parents=True, exist_ok=True)

        from censo_ext.Tools.utility import save_dict_orcaS
        save_dict_orcaS(Av_orcaS, self.Average_orcaSJ.orcaSParams)
        import json
        with open(Av_orcaAtoms, 'w') as f:
            f.write(json.dumps(self.Average_orcaSJ.idxAtoms))
        np.savetxt(Av_orcaJ, self.Average_orcaSJ.orcaJCoups, fmt="%10.5f")

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
        # print(int(nNuclei/6),nNuclei%6)
        nLines = 0
        for idx in range(int(nNuclei/6)):
            nLines += nNuclei - idx * 6 + 3
        end_idx = start_idx + nLines + nNuclei % 6 + 3 - 1
        # print(start_idx,end_idx)

        for idx, line in enumerate(lines):
            if idx >= start_idx and idx <= end_idx:
                DataJ.append(line.rstrip())
        # print("="*80)

        for idx, line in enumerate(lines):
            if re.search(r"\+\/\-", line):
                tmp: list[float] = [int(i) for i in line.split()[0:3]]
                tmp.append(float(line.split()[3]))
                self.anmrS.append(tmp)
        # print(self.anmrS)

        for idx, line in enumerate(lines):
            if re.search(r"1H resonance frequency", line):
                self.frq = float(line.split()[6])
        # print("self.anmrS",len(self.anmrS))
        # print(self.anmrS)
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
        # print(ListDataJDeleteBlank)

        # Delete the right serial numbers of JCoup
        ListDataDelete_nAtoms: list[str] = []
        for line in ListDataJDeleteBlank:
            if (len(ListDataDelete_nAtoms) % nLinesDataJ == 0):
                for j in range(6*((int)(len(ListDataDelete_nAtoms)/nLinesDataJ))):
                    ListDataDelete_nAtoms.append("\n")
            ListDataDelete_nAtoms.append(line[6:])
        # print(ListDataDelete_nAtoms)

        DataJ_triangle = [""]*(nLinesDataJ)

        # Restruct to one all JCoup table
        for idx, line in enumerate(ListDataDelete_nAtoms):
            DataJ_triangle[idx % nLinesDataJ] = DataJ_triangle[idx %
                                                               nLinesDataJ].rstrip("\n") + " " + line.rstrip("\n")
        # print(DataJ_triangle)
        self.anmrJ = np.zeros((nLinesDataJ, nLinesDataJ))

        # copy half to other half data on JCoup
        for idx in range(nLinesDataJ):
            for idy in range(idx):
                self.anmrJ[idy][idx] = self.anmrJ[idx][idy] = float(
                    DataJ_triangle[idx].split()[idy])
        # print(self.anmrJ)

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

    def method_read_nucinfo(self, filename: Path = Path("anmr_nucinfo")) -> None:
        filename = self.__Directory / Path(filename)
        is_exist(filename)

        lines: list[str] = open(filename, "r").readlines()
        nAtoms: int = int(lines[0].split()[0])
        # print("nAtoms : ", nAtoms)
        del lines[0]
        # print("page", len(lines)/nAtoms, nAtoms)
        page: list = []
        for idx, x in enumerate(lines):
            x = x.rstrip()
            # print(int((idx)/(nAtoms*2)),idx%2,idx,x)
            # page = int((idx)/(nAtoms*2))
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
            exit(0)

    def method_read_enso(self, fileName: Path = Path("anmr_enso")) -> None:
        ''' anmr_enso :  8 columns '''
        ''' ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi  '''

        fileName = self.__Directory / Path(fileName)
        is_exist(fileName)
        self.enso = np.genfromtxt(fileName, names=True)
        if len(self.enso[0]) != 8:
            print("something wrong in your anmr_enso file")
            ic()
            exit(0)

    def method_print_enso(self) -> None:
        print("ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi     ")
        for Enso in self.enso:
            print(f'{int(Enso[0]):<1d}      {int(Enso[1]):<4d} {int(Enso[2]):<4d} {Enso[3]:>6.4f} {Enso[4]: > 11.7f} {Enso[5]: > 10.7f} {Enso[6]: > 10.7f} {Enso[7]:>4.3f}')  # nopep8

    def method_save_enso(self, fileName: Path = Path("anmr_enso.new")) -> None:
        DirFileName: Path = self.__Directory / Path(fileName)
        original_stdout = sys.stdout
        with open(DirFileName, "w") as f:
            sys.stdout = f
            self.method_print_enso()
        sys.stdout = original_stdout


class OrcaSJ():
    def __init__(self) -> None:
        self.orcaJCoups: np.ndarray = np.array([])
        self.orcaSParams: dict[int, float] = {}
        self.orcaAnisotropy: dict[int, float] = {}
        self.CONFSerialNums: int = 0
        self.idxAtoms: dict[int, str] = {}

    def method_read_orcaJ(self, filename: Path = Path("orcaJ.out")) -> bool:
        ''' Read the file orcaJ.out in censo program '''
        print(" method_read_orcaJ", filename)
        is_exist(filename)

        start_idx, end_idx = 0, 0
        DataJ: list = []
        lines: list[str] = open(filename, "r").readlines()
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

        # end_idx = start_idx + nLines - 1
        # ic(start_idx,end_idx,nLines)
        if start_idx == 0 or end_idx == 0:
            print(filename, " the data of the file is some error ...")
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
        # ic(DataJ[0],len(DataJ[0]))
        self.orcaJCoups = np.array(DataJ, dtype=float)
        return True

    def method_read_orcaS(self, filename: Path = Path("orcaS.out")) -> bool:
        ''' Read the file orcaS.out in censo program '''
        print(" method_read_orcaS", filename)
        is_exist(filename)

        start_idx, end_idx = 0, 0
        DataS: list = []
        lines: list[str] = open(filename, "r").readlines()
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
            # ic(start_idx,end_idx,nNuclei)
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
        # ic(start_idx,end_idx)

        for idx, line in enumerate(lines):
            if idx >= start_idx and idx <= end_idx:
                DataS.append(line.rstrip())
        # ic(DataS)
        self.idxAtoms, self.orcaAnisotropy, self.orcaSParams = {}, {}, {}
        for idx, x in enumerate(DataS):
            # self.idxAtoms.append([int(x.split()[0])+1, x.split()[1]])
            self.idxAtoms[int(x.split()[0])+1] = str(x.split()[1])
            self.orcaSParams[int(x.split()[0])+1] = float(x.split()[2])
            self.orcaAnisotropy[int(x.split()[0])+1] = float(x.split()[3])
        # ic(self.orcaSParams)
        return True

    def method_save_orcaS(self) -> list:
        raise NotImplementedError
        if len(self.idxAtoms) == len(self.orcaSParams):
            lines: list = []
            lines.append("  Nucleus  Element    Isotropic     Anisotropy\n")
            lines.append("  -------  -------  ------------   ------------\n")

            for key, value in self.idxAtoms.items():
                str1: str = f"  \
                    {key-1: > 5d}{value: > 8s}{self.orcaSParams[key]: > 15.3f}{self.orcaAnisotropy[key]: > 15.3f}\n"
                lines.append(str1)
            return lines

        else:
            print("your orcaJ and orcaS is not fit each other")
            print("    exit and close the program !!! ")
            ic()
            exit(0)

    def method_save_orcaJ(self) -> list:
        raise NotImplementedError
        # print(" To save the orcaS.out and orcaJ.out from data")
        if len(self.idxAtoms) == len(self.orcaJCoups[0]):
            lines: list = []
            list_idxAtoms: list = list(self.idxAtoms.keys())

            for idx in range(len(self.idxAtoms)):
                for idy in range(idx+1, len(self.idxAtoms), 1):
                    # print(idx, idy)
                    str1: str = f" NUCLEUS A = {(self.idxAtoms[list_idxAtoms[idx]]): > s}\
                        {(list_idxAtoms[idx]-1): > 5d} NUCLEUS B = {self.idxAtoms[list_idxAtoms[idy]]: > s}{list_idxAtoms[idy]-1: > 5d}"
                    # " NUCLEUS A = H    4 NUCLEUS B = H    5"
                    str2: str = f" Total            0.000            0.000            0.000  iso = \
                        {self.orcaJCoups[idx][idy]: > 13.5f}"
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
        if len(self.idxAtoms) == len(self.orcaSParams):
            print("Nucleus  Element   Anisotropy")
            for idx, Atom in self.idxAtoms.items():
                print(f'{idx:>5d}', f'{Atom:>8s}', end="")
                print(f'{self.orcaSParams[idx]:>15.3f}')
        else:
            print("your orcaJ and orcaS is not fit each other")
            print("    exit and close the program !!! ")
            ic()
            exit(0)

    def method_print_orcaJ(self) -> None:
        for idx in range(self.orcaJCoups[0].size):
            for idy in range(self.orcaJCoups[0].size):
                print(f'{(self.orcaJCoups[idx][idy]):>8.3f}', end="")
            print("")


class CensoDat():

    def __init__(self, fileName: Path = Path("anmr.dat")) -> None:
        is_exist(fileName)
        self.__fileName: Path = Path(fileName)
        self.__dat: np.ndarray = np.genfromtxt(fileName)

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

        # ic(start, end, dpi)
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

    def set_fileName(self, fileName: Path) -> None:
        self.__fileName = Path(fileName)

    def get_fileName(self) -> Path:
        return self.__fileName

    def get_Dat(self) -> np.ndarray:
        return self.__dat
