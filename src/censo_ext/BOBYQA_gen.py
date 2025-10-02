#!/usr/bin/env python
from censo_ext.Tools.anmrfile import AD_BOBYQA, AD_Normal
from censo_ext.Tools.utility import cosine_similarity, print_arguments
import argparse
import numpy as np
import numpy.typing as npt
descr = """
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=f"{descr}",
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
    parser.add_argument(
        "-c",
        "--comb",
        dest="comb",
        action="store",
        required=False,
        type=int,
        default=16,
        help="Maximum of Combinations [default 16]",
    )

    parser.add_argument(
        "-l",
        "--limits",
        dest="limits",
        action="store",
        required=False,
        type=float,
        default=1,
        help="limits of border ppms [default 1]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    from censo_ext.Tools.datfile import Peaks_npz, unit_conversion
    blank: npt.NDArray = np.array([])
    uc = unit_conversion(blank)
    peaks = Peaks_npz(uc, args.file)
    print("  ========== Start ==========")
    peaks.method_read_file()
    peaks.method_print()

    AD_normal: AD_Normal = AD_Normal()

    if (AD_normal.Exist()):
        AD_normal.method_load_files()
        if isinstance(AD_normal.SParams, dict):
            OrcaS = np.array(list(AD_normal.SParams.items()))
            print(f"\n{AD_normal._file_orcaS}\n{OrcaS}")
            in_SParams: list[int] = list(map(int, AD_normal.SParams.keys()))
        else:
            print("The OrcaS.out in Averaage Directory is not dict")
            exit(1)
        from censo_ext.Tools.anmrfile import Anmr
        inAnmr: Anmr = Anmr()
        inAnmr.method_read_nucinfo()
        ChemEqvs: dict[int, list[int]] = {key: value for key, value in
                                          inAnmr.NeighborChemEqvs.items()
                                          if key in in_SParams}
        Groups: list[list[int]] = list(
            sorted(value) for value in ChemEqvs.values())
        unique_group: list[list[int]] = []
        for item in Groups:
            if item not in unique_group:
                unique_group.append(item)
        unique_group_first_idx: list[int] = []
        for item in unique_group:
            unique_group_first_idx.append(item[0])

        # Normal Sim_SParams/sorted_SParams is more than real known peaks
        Sim_SParams: list[float] = []
        for x in unique_group_first_idx:
            index = np.where(OrcaS.T[0] == x)
            Sim_SParams.append(-float(OrcaS.T[1][index][0]))

        sorted_SParams: npt.NDArray[np.float64] = np.array(sorted(Sim_SParams))
        real_peaks: npt.NDArray[np.float64] = peaks.get_cIDs_center_peaks()
        print(f"\nSorted SParams : \n{sorted_SParams}")
        # print(f"{len(sorted_SParams)=}")
        # print(f"{len(real_peaks.T)=}")
        if len(sorted_SParams) < len(real_peaks.T):
            print("  The numbers of real peaks are more than the numbers of calculation")
            exit(0)
        from censo_ext.Tools.utility import sub_numpy
        nGroups: npt.NDArray[np.int64] = sub_numpy(
            sorted_SParams, args.comb)
        OrcaS_BOBYQA = np.insert(OrcaS, 2, 0, axis=1)
        # ic(nGroups)
        c_intp = np.argmax(nGroups)
        c_start = np.sum(nGroups[:c_intp])
        c_end = np.sum(nGroups[:c_intp+1])-1
        c_start_ppm = sorted_SParams[c_start]
        c_end_ppm = sorted_SParams[c_end]
        # ic(c_start_ppm, c_end_ppm)

        ordered_orcaS: npt.NDArray[np.float64] = np.array([])
        counter: int = 1
        for idx, groups in enumerate(nGroups):

            if groups != np.max(nGroups):
                # ic(idx, groups)
                start = np.sum(nGroups[:idx])
                end = np.sum(nGroups[:idx+1])-1
                # ic(start, end)
                Wait_Check_SParams = sorted_SParams[start:end+1]
                if sorted_SParams[start] > c_start_ppm:
                    S_start = sorted_SParams[start]
                else:
                    S_start = sorted_SParams[start]-args.limits
                if sorted_SParams[start] < c_end_ppm:
                    S_end = sorted_SParams[end]
                else:
                    S_end = sorted_SParams[end]+args.limits
                a = np.argwhere(real_peaks[1] > S_start)
                b = np.argwhere(real_peaks[1] < S_end)
                total_set = np.array(
                    list(set(a.T[0]).intersection(set(b.T[0]))))
                Wait_Check_Reals = np.array(real_peaks[1][total_set])
                # ic(Wait_Check_Reals)
                # ic(Wait_Check_SParams)
                if len(Wait_Check_Reals) >= len(Wait_Check_SParams):
                    Large = Wait_Check_Reals
                    Small = Wait_Check_SParams
                if len(Wait_Check_SParams) > len(Wait_Check_Reals):
                    Large = Wait_Check_SParams
                    Small = Wait_Check_Reals

                from itertools import combinations
                combs = list(combinations(
                    [*range(0, len(Large))], len(Small)))  # type: ignore # nopep8
                result: list[float] = []
                print(f"The numbers of Combinations : {len(combs)}")
                # ic(combs)
                for comb in combs:
                    result.append(cosine_similarity(
                        Large[list(comb)], Small))  # type: ignore # nopep8
                print("The result of cosine_similarity :")
                print(f"{np.array(result)}")
                print(f"The best of the result : {max(result)}")
                args_combs: tuple[int, ...] = combs[result.index(max(result))]
                # ic(Large, args_combs) #nopep8
                if len(Wait_Check_Reals) >= len(Wait_Check_SParams):
                    ordered_orcaS = Large[list(args_combs)]  # type: ignore # nopep8
                if len(Wait_Check_SParams) > len(Wait_Check_Reals):
                    ordered_orcaS = Large[list(args_combs)]  # type: ignore # nopep8
                print(f"The best of the array : {ordered_orcaS}")

                for idx, x in enumerate(ordered_orcaS):
                    intp = np.where(OrcaS.T[1] == -x)
                    OrcaS_BOBYQA.T[1][intp] = -Wait_Check_Reals[idx]
                    OrcaS_BOBYQA.T[2][intp] = counter
                    counter += 1

        AD_bobyqa: AD_BOBYQA = AD_BOBYQA()
        AD_bobyqa.JCoups = AD_normal.JCoups
        AD_bobyqa.idx1Atoms = AD_normal.idx1Atoms
        AD_bobyqa.SParams = OrcaS_BOBYQA
        AD_bobyqa.method_save_files()
        print("\nThe data is saved to orcaS-BOBYQA.out file")
        print(f"{AD_bobyqa._file_orcaS}\n{AD_bobyqa.SParams}")
        print("  ========== End ==========")

    else:
        exit(0)


if __name__ == "__main__":
    main()
