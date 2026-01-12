#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
from pathlib import Path
from censo_ext.Tools.utility import AtomID, delete_all_files
from censo_ext.Tools.utility import IsExist_bool, print_arguments

descr = """
________________________________________________________________________________
| Purpose : some JCoup constant in Eqv. atom in anmr program are not averge values,
|           so overwrite the orcaJ.out file to get average JCoup constant
| Default : Overwrite the orcaJ.out and backup the old data to orcaJ.out.backup
| Recover : -r Copy the orcaJ.out.backup to orcaJ.out [default False]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=f"{descr}",
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


def Atom_Equivalent(file: Path | str = Path("anmrh.out")) -> list[list[AtomID]]:
    file = Path(file)
    from censo_ext.Tools.anmrfile import Anmr
    inAnmr: Anmr = Anmr()
    inAnmr.method_read_anmrSJ(file)
    SParams: list = inAnmr.anmrS
    inAnmr.method_read_nucinfo()
    AtomEqv: list[list[AtomID]] = []
    for x in [a[1] for a in SParams]:
        AtomEqv.append(inAnmr.NeighborMangetEqvs[x])
    return AtomEqv


def function_read_orcaJ(file: Path = Path("orcaJ.out")) -> npt.NDArray[np.float64]:
    from censo_ext.Tools.anmrfile import OrcaSJ
    single_orcaSJ = OrcaSJ()
    if single_orcaSJ.method_read_orcaJ(file):
        return single_orcaSJ.JCoups
    else:
        print("  orcaJ.out file have something wrong !!!")
        print("  Exit and Close the program !!!")
        exit(0)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

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
        idx_h_lines: list[int] = []
        for idx0, line in enumerate(lines):
            if r"h" in line:
                idx_h_lines.append(idx0)
            # if r"h" in line:
            #    idx_h_lines.append(idx0)
        idx0_h_lines: npt.NDArray[np.int64] = np.array(idx_h_lines) - 1
        np.set_printoptions(formatter={'float': '{:12.5f}'.format})

        AtomIDs_Eqv: list[list[AtomID]] = Atom_Equivalent("anmrh.out")

        for idx, x in enumerate(AtomIDs_Eqv):
            for idy, y in enumerate(x):
                AtomIDs_Eqv[idx][idy] = (idx0_h_lines + 1).tolist().index(y)

        for dirName in (dirNames):
            fileBackup: Path = Path(f"{dirName}/NMR/orcaJ.out.backup")
            orcaJfile: Path = Path(f"{dirName}/NMR/orcaJ.out")
            JCoups: npt.NDArray[np.float64]

            if IsExist_bool(fileBackup):
                JCoups = function_read_orcaJ(fileBackup)
            else:
                JCoups = function_read_orcaJ(orcaJfile)
                import shutil
                shutil.copyfile(orcaJfile, fileBackup)

            for i in range(len(AtomIDs_Eqv)-1, -1, -1):

                if (len(AtomIDs_Eqv[i]) != 1):
                    print(f"{AtomIDs_Eqv[i]}=")
                    JCoup_temp = np.mean(JCoups[(AtomIDs_Eqv[i])], axis=0)
                    JCoups[AtomIDs_Eqv[i]] = JCoup_temp
                    JCoups.transpose()[AtomIDs_Eqv[i]] = JCoup_temp

            for i in range(len(AtomIDs_Eqv)-1, -1, -1):
                if (len(AtomIDs_Eqv[i]) > 2):
                    for j in range(len(AtomIDs_Eqv[i])-1, -1, -1):
                        for k in range(len(AtomIDs_Eqv[i])-1, -1, -1):
                            JCoups[AtomIDs_Eqv[i][j], AtomIDs_Eqv[i][k]] = 0
            np.set_printoptions(formatter={'float': '{:12.5f}'.format})

            orcaJ_File = (dirName + '/NMR/orcaJ.out')
            delete_all_files(orcaJ_File)
            with open(orcaJ_File, 'w') as outfile:
                for i in range(0, len(idx0_h_lines)):
                    for j in range(i+1, len(idx0_h_lines)):
                        print(
                            f" NUCLEUS A = H {int(idx0_h_lines[i])} NUCLEUS B = H {int(idx0_h_lines[j])}", file=outfile)
                        print(
                            f" Total            0.000            0.000            0.000  iso= {str(JCoups[i][j]):.5f}", file=outfile)

            print(f" Directory of saved file: {orcaJ_File}")


if __name__ == "__main__":
    main()
