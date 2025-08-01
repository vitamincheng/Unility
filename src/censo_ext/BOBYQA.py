#!/usr/bin/env python3
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
|                                          [07.22.2025] vitamin.cheng@gmail.com
| BOBYQA.py
| Usages  : BOBYQA.py [options]
| [options]
| Dir      : -d the directory [default .]
| Ref      : -r the actual reference file [default 1r.dat]
| Limit    : -l limit border(ppm) [defalut 0.20]
| Package  : Tools
| Module   : ??.py
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=argparse.RawDescriptionHelpFormatter,
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
        help="Use external anmr execute file",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


FileName_OrcaS: Path = Path("Average/NMR/orcaS.out")
FileName_BOBYQA: Path = Path("Average/NMR/orcaS-BOBYQA.out")


def rosenbrock(x0) -> float:
    orcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileName_BOBYQA)

    if len(x0) == 1:
        # single peak
        for idx_key in idx_keys:
            orcaS_Table[idx_key][1] = x0[0]
    else:
        # group peaks
        for idx, idx_key in enumerate(idx_keys):
            orcaS_Table[idx_key][1] = x0[idx]

    np.savetxt(Directory/FileName_BOBYQA, orcaS_Table, fmt="%10d %10.5f %10d")
    orcaS_Table = np.delete(orcaS_Table, 2, axis=1)
    np.savetxt(Directory/FileName_OrcaS, orcaS_Table, fmt="%10d %10.5f")

    from censo_ext.Tools.anmrfile import CensoDat

    if prog is True:

        cwd: Path = Path(os.getcwd())
        os.chdir(Directory)
        import sys
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
        import censo_ext.anmr as anmr
        x: dict = {'out': 'output.dat', "dir": Directory, "json": None, 'mf': 500.0, 'lw': None, 'ascal': None, 'bscal': None, 'thr': None, 'thrab': 0.025,
                   'tb': 4, 'cutoff': 0.001, 'start': None, 'end': None, 'show': False, 'mss': 9, 'auto': True, 'average': True, 'bobyqa': False}
        import sys
        sys.stdout = open(os.devnull, 'w')
        anmr.main(args=argparse.Namespace(**x))
        sys.stdout = sys.__stdout__

        Dat_Cal: CensoDat = CensoDat(file=Directory/Path(x["out"]))
    else:
        print("Something wrong in your argument ")
        ic()
        raise ValueError("Something wrong in your argument")

    Dat_Ref: CensoDat = CensoDat(file=Path(Directory/Dat_fileName))
    Dat_Cal.method_normalize_dat()
    Dat_Ref.method_normalize_dat()
    Diff: CensoDat = Dat_Cal.method_subtract_dat(Dat_Ref)

    res = np.sum(np.square(Diff.get_Dat()))
    # ic(Res)
    return res


def Scan_single_Peak() -> None:
    import pybobyqa
    OrcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileName_BOBYQA)
    in_set: set[int] = set(OrcaS_Table.T[2].astype(int).tolist())
    in_set = {x for x in in_set if x <= 100 and x >= 1}
    # in_set.discard(0)
    for nSerial in in_set:
        ic(nSerial)
        intp: npt.NDArray[np.int64] = np.argwhere(
            OrcaS_Table.T[2] == nSerial).flatten()
        Data_Chemical_Shift: list[float] = list(
            map(float, np.atleast_1d(OrcaS_Table.T[1][intp[0]])))
        ic(Data_Chemical_Shift)
        ic(intp)

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
    orcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileName_BOBYQA)

    in_set: set[int] = set(orcaS_Table.T[2].astype(int).tolist())
    in_set = {x for x in in_set if x >= 1000}

    for nSerial in in_set:

        global idx_keys
        Data_Chemical_Shift: list[float] = list()

        ic(nSerial)
        intp: npt.NDArray[np.int64] = np.argwhere(
            orcaS_Table.T[2] == nSerial).flatten()
        Data_Chemical_Shift: list[float] = list(
            map(float, orcaS_Table.T[1][intp]))
        idx_keys = list(map(int, intp))
        idx_atoms = orcaS_Table.T[0][intp]
        ic(Data_Chemical_Shift)
        ic(intp)
        ic(idx_atoms)

        nNumbers: int = len(idx_keys)
        from itertools import permutations
        Permutations: list[tuple] = list(permutations(
            [*range(0, nNumbers)], nNumbers))
        solution_f: list = []
        solution_x0: list = []

        for Permutation in Permutations:
            x0: npt.NDArray[np.int64] = np.array(Data_Chemical_Shift)[
                list(Permutation)]
            # idx_atoms =[ int(x) for x in np.array(idx_atoms)[list(Permutation)]]

            lower = x0 - limit_border
            upper = x0 + limit_border
            # ic(idx_atoms, x0.tolist())
            soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
                lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
            print(soln)
            solution_f.append(soln.f)
            solution_x0.append(soln.x)

        ic(solution_f)
        argsmin: int = min(range(len(solution_f)), key=solution_f.__getitem__)
        list_x0: list[int] = solution_x0[argsmin]

        OrcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
            Directory/FileName_BOBYQA)

        for idx, idx_key in enumerate(idx_keys):
            OrcaS_Table[idx_key][1] = list_x0[idx]

        np.savetxt(Directory/FileName_BOBYQA,
                   OrcaS_Table, fmt="%10d %10.5f %10d")
        OrcaS_Table = np.delete(OrcaS_Table, 2, axis=1)
        np.savetxt(Directory/FileName_OrcaS, OrcaS_Table, fmt="%10d %10.5f")
    print(" ==== Finished group_peaks ====")


def Create_BOBYQA() -> tuple[bool, bool]:
    orcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
        Directory / FileName_OrcaS)
    orcaS_Table = np.insert(orcaS_Table, 2, 0, axis=1)
    np.savetxt(Directory / FileName_BOBYQA,
               orcaS_Table, fmt="%10d %10.5f %10d")
    print(" Create the orcaS-BOBYQA.out file ")
    print(" three column :         0 - Do nothing ")
    print("                     1~99 - Use BOBYQA single point to find the peak (use First Chemical Shift)")
    print("               Above 1000 - Use BOBYQA groups to find the peaks ")
    print(" Run this program again")
    return (False, True)


def main(args: argparse.Namespace = argparse.Namespace()) -> tuple[bool, bool]:
    '''

    return (bool,bool)
           FileName_BOBYQA is Exist ??
           Create_BOBYQA() is Work ??
    '''
    global Directory
    global Dat_fileName
    global limit_border
    global prog
    if args == argparse.Namespace():
        args = cml("")
    if args.dir:                            # default .
        Directory = Path(args.dir)
        ic(Directory)
    if args.ref:                            # default 1r.dat
        Dat_fileName = args.ref
        ic(Dat_fileName)
    if args.limit:                          # default 0.20 ppm
        limit_border = args.limit
    if args.prog:
        prog = args.prog
    else:
        prog = False

    if IsExist_return_bool(Directory / FileName_OrcaS):                 # type: ignore # nopep8
        if IsExist_return_bool(Directory / FileName_BOBYQA):            # type: ignore # nopep8
            ic("BOBYQA is exist")
            if prog is True:
                # ic("External")
                cwd: Path = Path(os.getcwd())
                os.chdir(Directory)  # type: ignore
                print(" Need to build the new CONF* system")
                print(" And copy your CONF* to /Backup/CONF*")
                print(" And create a new CONF1 (copy from /Backup/CONF1)")
                print(" Modify from /Average/NMR/orcaS.out")
                Res = input("Are you Sure to Continue ?? (Y/N)")
                if Res == "Y":
                    Ref_TMS: float = 31.820
                    subprocess.call("mkdir backup", shell=True)
                    subprocess.call("mv CONF* backup", shell=True)
                    subprocess.call("cp -r backup/CONF1/ .", shell=True)
                    import sys
                    template_inp: str = "CONF1/NMR/orcaS.out"
                    original_stdout = sys.stdout
                    with open(template_inp, "w") as f:
                        sys.stdout = f
                        print("--------------------------------")
                        print("CHEMICAL SHIELDING SUMMARY (ppm)")
                        print("--------------------------------")
                        print("")
                        print("")
                        print("  Nucleus  Element    Isotropic     Anisotropy")
                        print("  -------  -------  ------------   ------------")
                    sys.stdout = original_stdout

                    orcaS_Table: npt.NDArray[np.float64] = np.genfromtxt(
                        FileName_OrcaS)  # type: ignore
                    orcaS_Table.T[1] = orcaS_Table.T[1]+Ref_TMS
                    orcaS_Table.T[0] = orcaS_Table.T[0]-1
                    np.savetxt(Path("CONF1/NMR/orcaS-main.out"), orcaS_Table, fmt="%7d       H    %10.5f          0")  # type: ignore # nopep8

                    subprocess.call(
                        "cat CONF1/NMR/orcaS-main.out >> CONF1/NMR/orcaS.out", shell=True)
                    subprocess.call("rm CONF1/NMR/orcaS-main.out", shell=True)
                os.chdir(cwd)
            Scan_single_Peak()
            Scan_group_Peaks()
            if prog is True:
                cwd: Path = Path(os.getcwd())
                os.chdir(Directory)  # type: ignore
                # print(" Recover the data from backup directory")
                # Res = input("Are you Sure to Continue ?? (Y/N)")
                # if Res == "Y":
                subprocess.call("rm -rf CONF1", shell=True)
                subprocess.call("mv backup/CONF* .", shell=True)
                os.chdir(cwd)
            return (True, False)
        else:
            return Create_BOBYQA()  # (False,True)
    else:
        return (False, False)


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
