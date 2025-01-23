#!/usr/bin/env python3
from icecream import ic
import argparse
import os
import numpy as np
from pathlib import Path

fileName: Path = Path("Average/NMR/orcaS.out")
fileNameB: Path = Path("Average/NMR/orcaS-BOBYQA.out")


def rosenbrock(x0) -> float:
    orcaS_Table: np.ndarray = np.genfromtxt(fileNameB)

    for idx, idx_k in enumerate(idx_key):
        orcaS_Table[idx_k][1] = x0[idx]

    np.savetxt(fileNameB, orcaS_Table, fmt="%10d %10.5f %10d")

    import anmr
    x: dict = {'out': 'output.dat', 'mf': 500.0, 'lw': None, 'ascal': None, 'bscal': None, 'thr': None, 'thrab': 0.025,
               'tb': 4, 'cutoff': 0.001, 'start': None, 'end': None, 'show': False, 'mss': 9, 'auto': True, 'average': True, 'bobyqa': False}
    import sys
    sys.stdout = open(os.devnull, 'w')
    np_dat: np.ndarray = anmr.main(args=argparse.Namespace(**x))
    sys.stdout = sys.__stdout__

    from Tools.anmrfile import CensoDat
    Dat_Cal: CensoDat = CensoDat(fileName=x["out"])
    Dat_Ref: CensoDat = CensoDat(fileName=Path("1r.dat"))
    Dat_Cal.method_normalize_dat()
    Dat_Ref.method_normalize_dat()
    Diff: CensoDat = Dat_Cal.method_subtract_dat(Dat_Ref)
    Diff.set_fileName(Path("diff.dat"))
    Diff.method_save_dat()

    rsquare = np.sum(np.square(Diff.get_Dat()))
    return rsquare


def single_scan() -> None:
    import pybobyqa
    orcaS_Table: np.ndarray = np.genfromtxt(fileNameB)

    for idx, SParam in enumerate(orcaS_Table):

        if SParam[2] == 1:

            data: list = []
            data.append(SParam[1])
            global idx_key
            idx_key = []
            idx_key.append(idx)
            x0 = np.array(data)
            lower = x0 - 0.2
            upper = x0 + 0.2
            ic(int(SParam[0]), x0.tolist())
            soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
                lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
            print(soln)


def group_scan(inGroup: list[int] = []) -> None:
    import pybobyqa
    orcaS_Table: np.ndarray = np.genfromtxt(fileNameB)

    group: set = set()
    for idx, SParam in enumerate(orcaS_Table):
        if SParam[2] >= 100:
            group.add(SParam[2])

    list_group: list = list(group)
    list_group.sort()

    for idg in list_group:

        global idx_key
        idx_key = list()
        data: list[float] = list()
        idx_atoms: list[int] = list()

        for idx, SParam in enumerate(orcaS_Table):
            if SParam[2] == idg:
                idx_key.append(idx)
                idx_atoms.append(int(SParam[0]))
                data.append(SParam[1])

        nNumbers: int = len(idx_key)
        from itertools import permutations
        Permutations: list = list(permutations(
            [*range(0, nNumbers)], nNumbers))
        solution_f: list = []
        solution_x0: list = []

        for Permutation in Permutations:
            x0 = np.array(data)[list(Permutation)]
            # idx_atoms =[ int(x) for x in np.array(idx_atoms)[list(Permutation)]]

            lower = x0 - 0.2
            upper = x0 + 0.2
            ic(idx_atoms, x0.tolist())
            soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=(
                lower, upper), scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
            print(soln)
            solution_f.append(soln.f)
            solution_x0.append(soln.x)

        ic(solution_f)
        argsmin: int = min(range(len(solution_f)), key=solution_f.__getitem__)
        list_x0: list[int] = solution_x0[argsmin]

        orcaS_Table = np.genfromtxt(fileNameB)

        for idx, idx_k in enumerate(idx_key):
            orcaS_Table[idx_k][1] = list_x0[idx]

        np.savetxt(fileNameB, orcaS_Table, fmt="%10d %10.5f %10d")


def new() -> None:
    from Tools.utility import IsExistReturnBool
    if IsExistReturnBool(fileNameB):
        return
    else:
        orcaS_Table: np.ndarray = np.genfromtxt(fileName)
        orcaS_Table = np.insert(orcaS_Table, 2, 0, axis=1)
        ic(orcaS_Table)
        np.savetxt(fileNameB, orcaS_Table, fmt="%10d %10.5f %10d")
        print(" Create the orcaS-BOBYQA.out file ")
        print(" three column :        0 - Do nothing ")
        print("                       1 - Use BOBYQA single point to find the peak ")
        print("               Above 100 - Use BOBYQA groups to find the peaks ")
        print(" Run this program again")
        ic()
        os._exit(0)


def main() -> None:
    print("="*80)
    new()
    single_scan()
    group_scan()


if __name__ == "__main__":
    main()
