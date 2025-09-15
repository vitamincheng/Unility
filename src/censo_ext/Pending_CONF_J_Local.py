#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
# import matplotlib.pyplot as plt
from sys import argv as sysargv
# from icecream import ic
from pathlib import Path
from censo_ext.Tools.utility import delete_all_files
from censo_ext.Tools.utility import IsExist_bool

descr = """
________________________________________________________________________________
|                                          [08.17.2024] vitamin.cheng@gmail.com 
| Purpose : some JCoup constant is not average in Eqv. atom in anmr program    
|           So overwrite the orcaJ.out file for average JCoup constant
| Default : Overwrite orcaJ.out and backup the old data to orcaJ.out.backup           
| Recover : -r Copy the orcaJ.out.backup to orcaJ.out [default False] 
| Needed  : orcaJ.out in each CONF folder, coord, anmrh.out, anmr_nucinfo   
| Package : Tools                  
| Module  : anmrfile.py 
|______________________________________________________________________________
"""


def cml():
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-r",
        "--recover",
        dest="recover",
        action="store_true",
        help="COPY the orcaJ.out.backup to orca.out [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def Atom_Equivalent(file: Path | str = Path("anmrh.out")) -> list:
    file = Path(file)
    from censo_ext.Tools.anmrfile import Anmr
    anmr = Anmr()
    anmr.method_read_anmrSJ(file)
    DataJ: list = anmr.anmrS
    anmr.method_read_nucinfo()
    AtomEqv: list = []
    for x in [a[1] for a in DataJ]:
        AtomEqv.append(anmr.NeighborMangetEqvs[x])
    return AtomEqv


def function_read_orcaJ(file: Path = Path("orcaJ.out")) -> npt.NDArray[np.float64]:
    from censo_ext.Tools.anmrfile import OrcaSJ
    single_orcaSJ = OrcaSJ()
    if single_orcaSJ.method_read_orcaJ(file):
        return single_orcaSJ.JCoups
    else:
        raise ValueError("  Someting wrong in your orcaJ.out file")


if __name__ == "__main__":

    args: argparse.Namespace = cml()
    print(f"    provided arguments: {" ".join(sysargv)}")
    print(descr)

    if args.recover:
        Dir: Path = Path.cwd()
        print(f"Files and directories in {Dir} : ")
        dirNames: list[str] = [x for x in Dir.walk()][0][1]

        idx: int = 0
        while (idx != len(dirNames)):
            if (dirNames[idx].find('CONF') == -1):
                del dirNames[idx]
            else:
                idx += 1
        print(f" Directories = {dirNames}")

        for dirName in dirNames:
            PathBackup: str = f"{dirName}/NMR/orcaJ.out.backup"
            orcaJPath: str = f"{dirName}/NMR/orcaJ.out"

            if IsExist_bool(PathBackup):
                import shutil
                shutil.copyfile(PathBackup, orcaJPath)
                delete_all_files(PathBackup)
            else:
                raise FileNotFoundError(
                    " Something wrong in your orcaJ.out folder")

        print("  Copy orcaJ.out.backup to orcaJ.out in every NMR folder")
        print("  Recover the orcaJ.out file in your CONF folder")
        print("  Exit and Close the program !!!")
        exit(0)

    else:
        Dir: Path = Path.cwd()
        print(f"Files and directories in {Dir} : ")
        dirNames: list[str] = [x for x in Dir.walk()][0][1]

        idx: int = 0
        while (idx != len(dirNames)):
            if (dirNames[idx].find('CONF') == -1):
                del dirNames[idx]
            else:
                idx += 1
        print(f" Directories = {dirNames}")

        file: Path = Path("coord")
        file_exists: bool = IsExist_bool(file)
        if not file_exists:
            raise FileNotFoundError(f"{file} the file is not exist ...")

        lines: list[str] = open(file, "r").readlines()
        import re
        idx_h_lines: list[int] = []
        for idx0, line in enumerate(lines):
            if re.search(r"h", line):
                idx_h_lines.append(idx0)
        np_idx_h_lines: npt.NDArray[np.int64] = np.array(idx_h_lines) - 1
        np.set_printoptions(formatter={'float': '{:12.5f}'.format})

        idx_Atom_Eqv: list = Atom_Equivalent("anmrh.out")

        for idx, x in enumerate(idx_Atom_Eqv):
            for idy, y in enumerate(x):
                idx_Atom_Eqv[idx][idy] = (np_idx_h_lines + 1).tolist().index(y)

        for dirName in (dirNames):
            fileBackup: Path = Path(f"{dirName}/NMR/orcaJ.out.backup")
            orcaJfile: Path = Path(f"{dirName}/NMR/orcaJ.out")
            JCoup: npt.NDArray[np.float64]

            if IsExist_bool(fileBackup):
                JCoup = function_read_orcaJ(fileBackup)
            else:
                JCoup = function_read_orcaJ(orcaJfile)
                import shutil
                shutil.copyfile(orcaJfile, fileBackup)

            for i in range(len(idx_Atom_Eqv)-1, -1, -1):

                if (len(idx_Atom_Eqv[i]) != 1):
                    print(f"{idx_Atom_Eqv[i]}=")
                    JCoup_temp = np.mean(JCoup[(idx_Atom_Eqv[i])], axis=0)
                    JCoup[idx_Atom_Eqv[i]] = JCoup_temp
                    JCoup.transpose()[idx_Atom_Eqv[i]] = JCoup_temp

            for i in range(len(idx_Atom_Eqv)-1, -1, -1):
                if (len(idx_Atom_Eqv[i]) > 2):
                    for j in range(len(idx_Atom_Eqv[i])-1, -1, -1):
                        for k in range(len(idx_Atom_Eqv[i])-1, -1, -1):
                            JCoup[idx_Atom_Eqv[i][j], idx_Atom_Eqv[i][k]] = 0
            np.set_printoptions(formatter={'float': '{:12.5f}'.format})

            orcaJ_File = (dirName + '/NMR/orcaJ.out')
            delete_all_files(orcaJ_File)
            with open(orcaJ_File, 'w') as outfile:
                for i in range(0, len(np_idx_h_lines)):
                    for j in range(i+1, len(np_idx_h_lines)):
                        outfile.write(
                            f" NUCLEUS A = H {int(np_idx_h_lines[i])} NUCLEUS B = H {int(np_idx_h_lines[j])}\n")
                        outfile.write(
                            f" Total            0.000            0.000            0.000  iso= {str(JCoup[i][j]):.5f}\n")

            print(f" Directory of saved file: {dirName}/NMR/orcaJ.out")
