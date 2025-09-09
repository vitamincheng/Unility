#!/usr/bin/env python
from icecream import ic
import argparse
import os
import numpy as np
import numpy.typing as npt
from pathlib import Path
import subprocess
from censo_ext.Tools.utility import IsExist_return_bool

descr = """
________________________________________________________________________________
| BOBYQA.py
| Usages   : BOBYQA.py [options]
| [options]
| Dir      : -d the directory [default .]
| Ref      : -r the actual reference file [default 1r.dat]
| Limit    : -l limit border(ppm) [defalut 0.20]
| Prog     : -p --prog Use external anmr execute file [default False]
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-d",
        "--dir",
        dest="dir",
        action="store",
        required=False,
        default=".",
        help="Provide output_file name [default .]",
    )

    parser.add_argument(
        "-r",
        "--ref",
        dest="ref",
        action="store",
        required=False,
        default="1r.dat",
        help="Provide ref file(dat) name [1r.dat]",
    )

    parser.add_argument(
        "-l",
        "--limit",
        dest="limit",
        action="store",
        type=float,
        required=False,
        default=0.20,
        help="Provide limit border (ppm) [0.20]",
    )

    parser.add_argument(
        "-p",
        "--prog",
        dest="prog",
        action="store_true",
        help="Use external anmr execute file [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


FileOrcaS: Path = Path("Average/NMR/orcaS.out")
FileBOBYQA: Path = Path("Average/NMR/orcaS-BOBYQA.out")

Ref_TMS: float = 31.820


def rosenbrock(x0) -> float:
    #
    # Average/NMR/orcaS-BOBYQA.out for setting
    # Average/NMR/orcaS.out        for anmr.py   (internal)
    # CONF1/NMR/orcaS.out          for anmr      (external)
    #
    from censo_ext.Tools.datfile import CensoDat
    orcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileBOBYQA)

    if len(x0) == 1:
        # single peak
        for idx_key in idx_keys:
            orcaS_Table[idx_key][1] = x0[0]
    else:
        # group peaks
        for idx, loop in enumerate(idx_keys):
            for idx_key in loop:
                orcaS_Table[idx_key][1] = x0[idx]

    np.savetxt(Directory/FileBOBYQA, orcaS_Table, fmt="%10d %10.5f %10d")
    orcaS_Table = np.delete(orcaS_Table, 2, axis=1)

    if prog:
        ic("External program: anmr")
        import sys
        cwd: Path = Path(os.getcwd())
        os.chdir(Directory)
        template_inp: str = "CONF1/NMR/orcaS.out"
        with open(template_inp, "w") as f:
            sys.stdout = f
            print("--------------------------------")
            print("CHEMICAL SHIELDING SUMMARY (ppm)")
            print("--------------------------------")
            print("")
            print("")
            print("  Nucleus  Element    Isotropic     Anisotropy")
            print("  -------  -------  ------------   ------------")
        sys.stdout = sys.__stdout__
        orcaS_Table.T[1] = orcaS_Table.T[1] + Ref_TMS
        orcaS_Table.T[0] = orcaS_Table.T[0]-1
        np.savetxt(Path("CONF1/NMR/orcaS-main.out"), orcaS_Table, fmt="%7d       H    %10.5f          0")  # type: ignore # nopep8

        subprocess.call(
            "cat CONF1/NMR/orcaS-main.out >> CONF1/NMR/orcaS.out", shell=True)
        subprocess.call("rm CONF1/NMR/orcaS-main.out", shell=True)

        sys.stdout = open(os.devnull, 'w')
        result = subprocess.call("anmr.sh", shell=True)
        sys.stdout = sys.__stdout__

        if result != 0:
            ic("Cal.=================", result)
            raise ValueError(
                " call anmr.sh process have something wrong !!!")

        os.chdir(cwd)
        Dat_Cal: CensoDat = CensoDat(file=Directory/Path("anmr.dat"))

    elif not prog:
        # ic("Internal")
        np.savetxt(Directory/FileOrcaS, orcaS_Table, fmt="%10d %10.5f")
        import censo_ext.anmr as anmr
        x: dict = {'out': 'output.dat', "dir": Directory, "json": None, 'mf': 500.0,
                   'lw': None, 'ascal': None, 'bscal': None, 'thr': None, 'thrab': 0.025,
                   'tb': 4, 'cutoff': 0.001, 'start': None, 'end': None, 'show': False,
                   'mss': 9, 'auto': True, 'average': True, 'bobyqa': False}
        import sys
        sys.stdout = open(os.devnull, 'w')
        anmr.main(args=argparse.Namespace(**x))
        sys.stdout = sys.__stdout__

        Dat_Cal: CensoDat = CensoDat(file=Directory/Path(x["out"]))
    else:
        print("Something wrong in your argument ")
        ic()
        raise ValueError("Something wrong in your argument")

    Dat_Cal.method_normalize_dat()
    global Dat_Ref
    try:
        Diff: CensoDat = Dat_Cal - Dat_Ref
    except NameError:
        Dat_Ref = CensoDat(file=Path(Directory/Dat_fileName))
        Dat_Ref.method_normalize_dat()
        Diff: CensoDat = Dat_Cal - Dat_Ref

    # Dat_Ref = CensoDat(file=Path(Directory/Dat_fileName))
    # Dat_Ref.method_normalize_dat()
    # Dat_Cal.method_normalize_dat()
    # Diff: CensoDat = Dat_Cal - Dat_Ref

    res = np.sum(np.square(Diff.get_Dat()))

    return res


def Scan_single_Peak() -> None:
    import pybobyqa
    OrcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileBOBYQA)
    in_set: set[int] = set(OrcaS_Table.T[2].astype(int).tolist())
    in_set = {x for x in in_set if x < 1000 and x >= 1}
    # in_set.discard(0)
    for nSerial in in_set:
        # ic(nSerial)
        intp: npt.NDArray[np.int64] = np.argwhere(
            OrcaS_Table.T[2] == nSerial).flatten()
        Data_Chemical_Shift: list[float] = list(
            map(float, np.atleast_1d(OrcaS_Table.T[1][intp[0]])))
        # ic(Data_Chemical_Shift)
        # ic(intp)

        global idx_keys
        idx_keys = list(intp)
        x0: npt.NDArray[np.float64] = np.array(Data_Chemical_Shift)
        lower: npt.NDArray[np.float64] = x0 - limit_border
        upper: npt.NDArray[np.float64] = x0 + limit_border
        # ic(int(SParam[0]), x0.tolist())
        soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
            lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
        print(soln)
    print(" ==== Finished single_peak ====")


def Scan_group_Peaks() -> None:
    import pybobyqa
    OrcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileBOBYQA)

    in_set: set[int] = set(OrcaS_Table.T[2].astype(int).tolist())
    in_set = {x for x in in_set if x >= 1000}

    if len(in_set) == 0:
        print(" ==== Finished group_peaks ====")
        return

    # Data structure of idx, Chemical_Shift, idx_atoms
    Data: list[list] = []

    for nSerial in in_set:
        intp: npt.NDArray[np.int64] = np.argwhere(
            OrcaS_Table.T[2] == nSerial).flatten()
        # idx_keys = list(map(int, intp))
        # idx_atoms = OrcaS_Table.T[0][intp]
        Data.append(
            [nSerial, OrcaS_Table.T[1][intp[0]], intp])

    # ic(Data)
    nNumbers: int = len(Data)
    from itertools import permutations
    Permutations: list[tuple] = list(permutations(
        [*range(0, nNumbers)], nNumbers))
    solution_f: list = []
    solution_x0: list = []
    global idx_keys

    for Permutation in Permutations:
        x0: npt.NDArray[np.float64] = np.array([x[1] for x in Data])[
            list(Permutation)]
        idx_keys = [x[2] for x in Data]
        # ic(idx_keys)
        lower = x0 - limit_border
        upper = x0 + limit_border
        soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
            lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
        print(soln)
        solution_f.append(soln.f)
        solution_x0.append(soln.x)

    ic(solution_f)
    argsmin: int = min(range(len(solution_f)), key=solution_f.__getitem__)
    list_x0: list[float] = solution_x0[argsmin]
    ic(list_x0)

    # After Permutations, the best choice x0 is calculated again and get output.dat or anmr.dat
    x0 = np.array(list_x0)
    limit_tiny: float = 0.0001
    lower = x0 - limit_tiny
    upper = x0 + limit_tiny
    soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
        lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)

    print(" ==== Finished group_peaks ====")


def Create_BOBYQA() -> None:
    OrcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileOrcaS)
    OrcaS_Table = np.insert(OrcaS_Table, 2, 0, axis=1)
    np.savetxt(Directory / FileBOBYQA,
               OrcaS_Table, fmt="%10d %10.5f %10d")
    print(" Create the orcaS-BOBYQA.out file ")
    print(" three column :         0 - Do nothing ")
    print("                     1~99 - Use BOBYQA to calcuate and fit each chemical shift of each number")
    print("               Above 1000 - Use BOBYQA to calucate and fit each chemical shift of all groups to find the peaks in one time")
    print("                            Each chemical shift of the same number is assigned to the same by the first chemical shift")
    print(" Run this program again")
    import sys
    sys.exit(0)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    '''

    return (bool,bool)
    '''
    #
    # orcaS-BOBYQA.out  for setting
    # Average/NMR/orcaS.out for anmr.py   (internal)
    # CONF1/NMR/orcaS.out for anmr        (external)
    #

    global Directory
    global Dat_fileName
    global limit_border
    global prog

    if args == argparse.Namespace():
        args = cml("")
    if args.dir:                            # default .
        Directory = Path(args.dir)
        # ic(Directory)
    if args.ref:                            # default 1r.dat
        Dat_fileName = args.ref
        # ic(Dat_fileName)
    if args.limit:                          # default 0.20 ppm
        limit_border = args.limit
    if args.prog:
        prog = args.prog
    else:
        prog = False

    if IsExist_return_bool(Directory / FileOrcaS):                 # type: ignore # nopep8
        ic(FileOrcaS, " is exist")
        if IsExist_return_bool(Directory / FileBOBYQA):            # type: ignore # nopep8
            ic(FileBOBYQA, " is exist")
            if prog:
                # ic("External")
                cwd: Path = Path(os.getcwd())
                os.chdir(Directory)  # type: ignore
                print(" Need to build the new CONF* system")
                print(" And copy your CONF* to /Backup/CONF*")
                print(" And create a new CONF1 (copy from /Backup/CONF1)")
                print(" Modify from /Average/NMR/orcaS.out")
                Res = input("Are you Sure to Continue ?? (Y/N)")
                if Res == "Y":
                    subprocess.call("mkdir backup", shell=True)
                    subprocess.call("mv CONF* backup", shell=True)
                    subprocess.call("cp -r backup/CONF1/ .", shell=True)
                os.chdir(cwd)
            Scan_single_Peak()
            Scan_group_Peaks()
            if prog:
                cwd: Path = Path(os.getcwd())
                os.chdir(Directory)  # type: ignore
                # print(" Recover the data from backup directory")
                # Res = input("Are you Sure to Continue ?? (Y/N)")
                # if Res == "Y":
                subprocess.call("rm -rf CONF1", shell=True)
                subprocess.call("mv backup/CONF* .", shell=True)
                subprocess.call("rmdir backup", shell=True)
                os.chdir(cwd)
        else:
            Create_BOBYQA()
    else:
        ic()
        raise FileNotFoundError(
            str(Directory/FileOrcaS) + " is not exist !!!")  # type: ignore # nopep8


if __name__ == "__main__":
    main()

# standard test
# BOBYQA.py -d tests/data/06.EthylAcetate/03.Censo -r 1r.dat
# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) -r 1r_h.dat
#
# Want to use external anmr program (Not Implemented on pytest,because subprocess function)
# Start new one
# Run below anmr.py command and it will create average folder
# anmr.py   -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) --auto
# Run below BOBYQA.py command and it will create ocaS-BOBYQA.out
# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS)
# rewrite orcaS-BOBYQA.out file and run below BOBYQA.py by use external program
# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) -r 1r_h.dat -p
