#!/usr/bin/env python
from censo_ext.Tools.utility import IsExist_bool, cosine_similarity
import argparse
import numpy as np
import numpy.typing as npt
from pathlib import Path
descr = """
________________________________________________________________________________
| For Generation of OrcaS.BOBYQA
| Usages   : BOBYQA_gen.py <geometry> [options]
| [options]
| input    : -i the input npz file [default peaks.npz]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="peaks.npz",
        help="Provide one input xyz file [default peaks.npz]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    # print_descr(descr)

    from censo_ext.Tools.datfile import Peaks_npz, unit_conversion
    blank: npt.NDArray = np.array([])
    uc = unit_conversion(blank)
    peaks = Peaks_npz(uc, args.file)
    print("  ========== Start ==========")
    peaks.method_read_file()
    peaks.method_print()

    FileOrcaS: Path = Path("Average/NMR/orcaS.out")
    if IsExist_bool(FileOrcaS):
        OrcaS: npt.NDArray[np.float64] = np.genfromtxt(FileOrcaS)
        print(f"\n{FileOrcaS}\n{OrcaS}")

        in_S: list[int] = list(map(int, OrcaS.T[0]))
        from censo_ext.Tools.anmrfile import Anmr
        inAnmr: Anmr = Anmr()
        inAnmr.method_read_nucinfo()
        ChemEqvs: dict[int, list[int]] = {key: value for key,
                                          value in inAnmr.NeighborChemEqvs.items() if key in in_S}
        Groups: list[list[int]] = list(
            sorted(value) for value in ChemEqvs.values())
        unique_group: list[list[int]] = []
        for item in Groups:
            if item not in unique_group:
                unique_group.append(item)
        unique_group_first_idx: list[int] = []
        for item in unique_group:
            unique_group_first_idx.append(item[0])

        # Normal Sim_SParams is more than real known peaks
        Sim_SParams: list[float] = []
        for x in unique_group_first_idx:
            index = np.where(OrcaS.T[0] == x)
            Sim_SParams.append(-float(OrcaS.T[1][index][0]))

        sorted_SParams: npt.NDArray[np.float64] = np.array(sorted(Sim_SParams))
        real_peaks: npt.NDArray[np.float64] = peaks.get_cIDs_center_peaks()
        print(sorted_SParams)

        from itertools import combinations
        combs = list(combinations(
            [*range(0, len(sorted_SParams))], len(real_peaks.T)))
        result: list[float] = []
        print(f"\nThe numbers of Combinations : {len(combs)}")
        for comb in combs:
            result.append(cosine_similarity(
                real_peaks[1], sorted_SParams[list(comb)]))
        print("The result of cosine_similarity :")
        print(f"{np.array(result)}")
        print(f"The best of the result : {max(result)}")
        args_combs: tuple[int, ...] = combs[result.index(max(result))]
        ordered_orcaS: npt.NDArray[np.float64] = sorted_SParams[list(
            args_combs)]
        print(f"The best of the array : {ordered_orcaS}")
        OrcaS_BOBYQA = np.insert(OrcaS, 2, 0, axis=1)
        for idx, x in enumerate(ordered_orcaS):
            intp = np.where(OrcaS.T[1] == -x)
            OrcaS_BOBYQA.T[1][intp] = -real_peaks[1][idx]
            OrcaS_BOBYQA.T[2][intp] = idx+1

        FileBOBYQA: Path = Path("Average/NMR/orcaS-BOBYQA.out")
        np.savetxt(FileBOBYQA, OrcaS_BOBYQA, fmt="%10d   %10.5f %10d")
        print("\nThe data is saved to orcaS-BOBYQA.out file")
        print(f"{FileBOBYQA}\n{OrcaS_BOBYQA}")
        print("  ========== End ==========")

    else:
        exit(0)


if __name__ == "__main__":
    main()
