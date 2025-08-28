#!/usr/bin/env python
import os
import argparse
import numpy as np
import numpy.typing as npt
# import matplotlib.pyplot as plt
from sys import argv as sysargv
from icecream import ic
from pathlib import Path

descr = """
________________________________________________________________________________
|                                          [08.17.2024] vitamin.cheng@gmail.com 
| Purpose : some JCoup constant is not average in Eqv. atom in anmr program    
|           So overwrite the orcaJ.out file for average JCoup constant
| Default : Overwrite orcaJ.out and backup the old data to orcaJ.out.backup           
| Recover : -r Copy the orcaJ.out.backup to orcaJ.out  
| Needed  : orcaJ.out in each CONF folder, coord, anmrh.out, anmr_nucinfo   
| Package : Tools                  
| Module  : anmrfile.py 
|______________________________________________________________________________
"""


def cml():
    # def cml(descr):
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        #        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )  # argparse.RawDescriptionHelpFormatter) #,

    parser.add_argument(
        "-r",
        "--recover",
        dest="recover",
        action="store_true",
        help="COPY the orcaJ.out.backup to orca.out ",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def Atom_Equivalent(file: Path = Path("anmrh.out")) -> list:
    from censo_ext.Tools.anmrfile import Anmr
    anmr = Anmr()
    anmr.method_read_anmrSJ(file)
    DataJ: list = anmr.anmrS
    anmr.method_read_nucinfo()
    AtomEqv: list = []
    for x in ([a[1] for a in DataJ]):
        # AtomEqv.append(anmr.nucinfo[1][x-1][2])
        AtomEqv.append(anmr.NeighborMangetEqvs[x])
    return AtomEqv


def function_read_orcaJ(file: Path = Path("orcaJ.out")) -> npt.NDArray[np.float64]:
    from censo_ext.Tools.anmrfile import OrcaSJ
    single_orcaSJ = OrcaSJ()
    if single_orcaSJ.method_read_orcaJ(file):
        return single_orcaSJ.JCoups
    else:
        print("  Someting wrong in your orcaJ.out file")
        print("  Exit and Close the program !!!")
        ic()
        exit(1)


if __name__ == "__main__":

    args: argparse.Namespace = cml()
    print("    provided arguments: {}".format(" ".join(sysargv)))
    print(descr)

    if args.recover:
        path = os.getcwd()
        print("Files and directories in ", path, " : ")
        dirNames: list[str] = next(os.walk(path))[1]

        idx: int = 0
        while (idx != len(dirNames)):
            if (dirNames[idx].find('CONF') == -1):
                del dirNames[idx]
            else:
                idx = idx+1
        print(" Directories = ", dirNames)

        for dirName in (dirNames):
            PathBackup: str = dirName + "/NMR/orcaJ.out.backup"
            orcaJPath: str = dirName + "/NMR/orcaJ.out"

            if os.path.exists(PathBackup):
                import shutil
                shutil.copyfile(PathBackup, orcaJPath)
                os.remove(PathBackup)
            else:
                print(" Something wrong in your orcaJ.out folder")
                print("", PathBackup)
                print("  Exit and Close the program !!!")
                ic()
                exit(1)

        print("  Copy orcaJ.out.backup to orcaJ.out in every NMR folder")
        print("  Recover the orcaJ.out file in your CONF folder")
        print("  Exit and Close the program !!!")
        exit(0)

    else:
        path: str = os.getcwd()
        print("Files and directories in ", path, " : ")
        dirNames: list[str] = next(os.walk(path))[1]

        idx = 0
        while (idx != len(dirNames)):
            if (dirNames[idx].find('CONF') == -1):
                del dirNames[idx]
            else:
                idx = idx+1
        print(" Directories = ", dirNames)

        from os.path import exists
        file: str = "coord"
        file_exists: bool = exists(file)
        if not file_exists:
            print(file, " the file is not exist ...")
            print("  Exit and Close the program !!!")
            ic()
            exit(1)

        lines: list[str] = open(file, "r").readlines()
        import re
        idx_h_lines: list[int] = []
        for idx, line in enumerate(lines):
            if re.search(r"h", line):
                idx_h_lines.append(idx)
        np_idx_h_lines: npt.NDArray[np.int64] = np.array(idx_h_lines) - 1
        np.set_printoptions(formatter={'float': '{:12.5f}'.format})

        idx_Atom_Eqv: list = Atom_Equivalent(Path("anmrh.out"))

        for idx, x in enumerate(idx_Atom_Eqv):
            for idy, y in enumerate(x):
                idx_Atom_Eqv[idx][idy] = (np_idx_h_lines + 1).tolist().index(y)

        for dirName in (dirNames):
            fileBackup: Path = Path(dirName + "/NMR/orcaJ.out.backup")
            orcaJfile: Path = Path(dirName + "/NMR/orcaJ.out")
            JCoup: npt.NDArray[np.float64]

            if os.path.exists(fileBackup):
                JCoup = function_read_orcaJ(fileBackup)
            else:
                JCoup = function_read_orcaJ(orcaJfile)
                import shutil
                shutil.copyfile(orcaJfile, fileBackup)

            for i in range(len(idx_Atom_Eqv)-1, -1, -1):

                if (len(idx_Atom_Eqv[i]) != 1):
                    ic(idx_Atom_Eqv[i])
                    JCoup_temp = np.mean(JCoup[(idx_Atom_Eqv[i])], axis=0)
                    JCoup[idx_Atom_Eqv[i]] = JCoup_temp
                    JCoup.transpose()[idx_Atom_Eqv[i]] = JCoup_temp

            for i in range(len(idx_Atom_Eqv)-1, -1, -1):
                if (len(idx_Atom_Eqv[i]) > 2):
                    for j in range(len(idx_Atom_Eqv[i])-1, -1, -1):
                        for k in range(len(idx_Atom_Eqv[i])-1, -1, -1):
                            JCoup[idx_Atom_Eqv[i][j], idx_Atom_Eqv[i][k]] = 0
            np.set_printoptions(formatter={'float': '{:12.5f}'.format})

            os.remove(dirName + '/NMR/orcaJ.out')
            with open(dirName + '/NMR/orcaJ.out', 'w') as outfile:
                for i in range(0, len(np_idx_h_lines)):
                    for j in range(i+1, len(np_idx_h_lines)):
                        outfile.write(
                            " NUCLEUS A = H "+str(int(np_idx_h_lines[i]))+" NUCLEUS B = H "+str(int(np_idx_h_lines[j]))+"\n")
                        outfile.write(" Total            0.000            0.000            0.000  iso= "+str(
                            format(JCoup[i][j], '.5f'))+"\n")

            print(" Directory of saved file: ", dirName + '/NMR/orcaJ.out')
