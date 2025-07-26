#!/usr/bin/env python3
from icecream import ic
import argparse
import os
import numpy as np
from pathlib import Path


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

    args: argparse.Namespace = parser.parse_args()
    return args


FileName_OrcaS: Path = Path("Average/NMR/orcaS.out")
FileName_BOBYQA: Path = Path("Average/NMR/orcaS-BOBYQA.out")


def rosenbrock(x0) -> float:
    orcaS_Table: np.ndarray = np.genfromtxt(Directory / FileName_BOBYQA)

    for idx, idx_k in enumerate(idx_keys):
        orcaS_Table[idx_k][1] = x0[idx]

    np.savetxt(Directory/FileName_BOBYQA, orcaS_Table, fmt="%10d %10.5f %10d")

    import censo_ext.anmr as anmr
    x: dict = {'out': 'output.dat', "dir": Directory, "json": None, 'mf': 500.0, 'lw': None, 'ascal': None, 'bscal': None, 'thr': None, 'thrab': 0.025,
               'tb': 4, 'cutoff': 0.001, 'start': None, 'end': None, 'show': False, 'mss': 9, 'auto': True, 'average': True, 'bobyqa': False}
    import sys
    sys.stdout = open(os.devnull, 'w')
    np_dat: np.ndarray = anmr.main(args=argparse.Namespace(**x))
    sys.stdout = sys.__stdout__

    from censo_ext.Tools.anmrfile import CensoDat
    Dat_Cal: CensoDat = CensoDat(file=Directory/Path(x["out"]))
    Dat_Ref: CensoDat = CensoDat(file=Path(Directory/Dat_fileName))
    Dat_Cal.method_normalize_dat()
    Dat_Ref.method_normalize_dat()
    Diff: CensoDat = Dat_Cal.method_subtract_dat(Dat_Ref)
    Diff.set_fileName(Path("diff.dat"))
    Diff.method_save_dat()

    rsquare = np.sum(np.square(Diff.get_Dat()))
    return rsquare


def Scan_single_Peak() -> None:
    import pybobyqa
    OrcaS_Table: np.ndarray = np.genfromtxt(Directory / FileName_BOBYQA)

    for idx, SParam in enumerate(OrcaS_Table):

        if SParam[2] == 1:  # on/off

            Data_Chemical_Shift: list = []
            Data_Chemical_Shift.append(SParam[1])
            global idx_keys
            idx_keys = []
            idx_keys.append(idx)
            x0: np.ndarray = np.array(Data_Chemical_Shift)
            lower: np.ndarray = x0 - limit_border
            upper: np.ndarray = x0 + limit_border
            ic(int(SParam[0]), x0.tolist())
            soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
                lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
            print(soln)


def Scan_group_Peaks(inGroup: list[int] = []) -> None:
    import pybobyqa
    orcaS_Table: np.ndarray = np.genfromtxt(Directory / FileName_BOBYQA)

    group: set = set()
    for idx, SParam in enumerate(orcaS_Table):
        if SParam[2] >= 100:
            group.add(SParam[2])

    list_group: list = list(group)
    list_group.sort()

    for idg in list_group:

        global idx_keys
        idx_keys = list()
        Data_Chemical_Shift: list[float] = list()
        idx_atoms: list[int] = list()

        for idx, SParam in enumerate(orcaS_Table):
            if SParam[2] == idg:
                idx_keys.append(idx)
                idx_atoms.append(int(SParam[0]))
                Data_Chemical_Shift.append(SParam[1])

        nNumbers: int = len(idx_keys)
        from itertools import permutations
        Permutations: list = list(permutations(
            [*range(0, nNumbers)], nNumbers))
        solution_f: list = []
        solution_x0: list = []

        for Permutation in Permutations:
            x0 = np.array(Data_Chemical_Shift)[list(Permutation)]
            # idx_atoms =[ int(x) for x in np.array(idx_atoms)[list(Permutation)]]

            lower = x0 - limit_border
            upper = x0 + limit_border
            ic(idx_atoms, x0.tolist())
            soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
                lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
            print(soln)
            solution_f.append(soln.f)
            solution_x0.append(soln.x)

        ic(solution_f)
        argsmin: int = min(range(len(solution_f)), key=solution_f.__getitem__)
        list_x0: list[int] = solution_x0[argsmin]

        OrcaS_Table = np.genfromtxt(Directory/FileName_BOBYQA)

        for idx, idx_k in enumerate(idx_keys):
            OrcaS_Table[idx_k][1] = list_x0[idx]

        np.savetxt(Directory/FileName_BOBYQA,
                   OrcaS_Table, fmt="%10d %10.5f %10d")


def Create_BOBYQA() -> None:
    orcaS_Table: np.ndarray = np.genfromtxt(Directory / FileName_OrcaS)
    orcaS_Table = np.insert(orcaS_Table, 2, 0, axis=1)
    ic(orcaS_Table)
    np.savetxt(Directory / FileName_BOBYQA,
               orcaS_Table, fmt="%10d %10.5f %10d")
    print(" Create the orcaS-BOBYQA.out file ")
    print(" three column :        0 - Do nothing ")
    print("                       1 - Use BOBYQA single point to find the peak ")
    print("               Above 100 - Use BOBYQA groups to find the peaks ")
    print(" Run this program again")


def main(args: argparse.Namespace = argparse.Namespace()):

    global Directory
    global Dat_fileName
    global limit_border
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

    from censo_ext.Tools.utility import IsExist_return_bool
    if IsExist_return_bool(Directory / FileName_OrcaS):                 # type: ignore # nopep8
        if IsExist_return_bool(Directory / FileName_BOBYQA):            # type: ignore # nopep8
            ic("BOBYQA is exist")
            Scan_single_Peak()
            Scan_group_Peaks()
        else:
            Create_BOBYQA()


if __name__ == "__main__":
    main()

# test
# BOBYQA.py -d tests/data/EthylAcetate/03.Censo -r 1r.dat
# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) -r 1r_h.dat
