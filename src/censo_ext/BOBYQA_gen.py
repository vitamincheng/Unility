#!/usr/bin/env python
from censo_ext.Tools.anmrfile import AD_BOBYQA, AD_Normal
from censo_ext.Tools.utility import cosine_similarity, print_arguments
import argparse
import numpy as np
import numpy.typing as npt
descr = """
_______________________________________________________________________________
| For generate orcaS-BOBYQA.out   
| 
|______________________________________________________________________________
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
    uc = unit_conversion(np.array([]))
    peaks = Peaks_npz(uc, args.file)
    print("  ========== Start ==========")
    peaks.method_read_file()
    peaks.method_print()

    AD_orcaS: AD_Normal = AD_Normal()

    if (AD_orcaS.Exist()):
        AD_orcaS.method_load_files()
        if isinstance(AD_orcaS.SParams, dict):
            # convert dict to npt.NDArray to OrcaS
            OrcaS: npt.NDArray[np.float64] = np.array(
                list(AD_orcaS.SParams.items()))
            print(f"\n{AD_orcaS._file_orcaS}\n{OrcaS}")
            in_SParams: list[int] = list(map(int, AD_orcaS.SParams.keys()))
        else:
            print("  The OrcaS.out in Averaage Directory is not dict format !!!")
            print("  Exit and Close the program !!!")
            exit(0)

        from censo_ext.Tools.anmrfile import Anmr
        inAnmr: Anmr = Anmr()
        inAnmr.method_read_nucinfo()
        ChemEqvs: dict[int, list[int]] = {key: value for key, value in
                                          inAnmr.NeighborChemEqvs.items()
                                          if key in in_SParams}
        Sorted_ChemEqvs: list[list[int]] = list(
            sorted(value) for value in ChemEqvs.values())

        unique_ChemEqvs: list[list[int]] = []

        for Sorted_ChemEqv in Sorted_ChemEqvs:
            if Sorted_ChemEqv not in unique_ChemEqvs:
                unique_ChemEqvs.append(Sorted_ChemEqv)

        unique_ChemEqvs_first_idx: list[int] = []
        for item in unique_ChemEqvs:
            unique_ChemEqvs_first_idx.append(item[0])

        # Normal Sim_SParams/sorted_SParams is more than real known peaks
        Sim_SParams: list[float] = []
        for x in unique_ChemEqvs_first_idx:
            index = np.where(OrcaS.T[0] == x)
            Sim_SParams.append(-float(OrcaS.T[1][index][0]))

        sorted_SParams: npt.NDArray[np.float64] = np.array(sorted(Sim_SParams))
        real_Peaks: npt.NDArray[np.float64] = peaks.get_cIDs_center_peaks()
        print(f"\nSorted SParams : \n{sorted_SParams}")
        # print(f"{len(sorted_SParams)=}")
        # print(f"{len(real_peaks.T)=}")
        if len(sorted_SParams) < len(real_Peaks.T):
            print(
                "  The numbers of real peaks are more than the numbers of simulation in normal")
            print("  Exit and Close the program !!!")
            exit(0)
        from censo_ext.Tools.utility import sub_numpy
        nGroups: npt.NDArray[np.int64] = sub_numpy(
            sorted_SParams, args.comb)
        OrcaS_BOBYQA: npt.NDArray[np.float64] = np.insert(OrcaS, 2, 0, axis=1)
        # ic(nGroups)
        c_intp: np.intp = np.argmax(nGroups)
        c_start_ppm: np.float64 = sorted_SParams[np.sum(nGroups[:c_intp])]
        c_end_ppm: np.float64 = sorted_SParams[np.sum(nGroups[:c_intp+1])-1]
        # ic(c_start_ppm, c_end_ppm)

        ordered_orcaS: npt.NDArray[np.float64] = np.array([])
        counter: int = 1
        for idx, groups in enumerate(nGroups):

            if groups != np.max(nGroups):
                # ic(idx, groups)
                start: np.int64 = np.sum(nGroups[:idx])
                end: np.int64 = np.sum(nGroups[:idx+1])-1
                # ic(start, end)
                Wait_Check_SParams: npt.NDArray[np.float64] = sorted_SParams[start:end+1]
                if sorted_SParams[start] > c_start_ppm:
                    S_start = sorted_SParams[start]
                else:
                    S_start = sorted_SParams[start]-args.limits
                if sorted_SParams[start] < c_end_ppm:
                    S_end = sorted_SParams[end]
                else:
                    S_end = sorted_SParams[end]+args.limits
                a: npt.NDArray[np.intp] = np.argwhere(real_Peaks[1] > S_start)
                b: npt.NDArray[np.intp] = np.argwhere(real_Peaks[1] < S_end)
                total_set: npt.NDArray[np.int64] = np.array(
                    list(set(a.T[0]).intersection(set(b.T[0]))))
                Wait_Check_Reals: npt.NDArray[np.float64] = np.array(
                    real_Peaks[1][total_set])
                # ic(Wait_Check_Reals)
                # ic(Wait_Check_SParams)
                if len(Wait_Check_Reals) >= len(Wait_Check_SParams):
                    Large = Wait_Check_Reals
                    Small = Wait_Check_SParams
                if len(Wait_Check_SParams) > len(Wait_Check_Reals):
                    Large = Wait_Check_SParams
                    Small = Wait_Check_Reals

                from itertools import combinations
                combs: list[tuple[int, ...]] = list(combinations(
                    [*range(0, len(Large))], len(Small)))  # type: ignore # nopep8
                result: list[float] = []
                print(f"The numbers of Combinations : {len(combs)}")
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
                    OrcaS_BOBYQA.T[1][intp] = - Wait_Check_Reals[idx]
                    OrcaS_BOBYQA.T[2][intp] = counter
                    counter += 1

        AD_bobyqa: AD_BOBYQA = AD_BOBYQA()
        AD_bobyqa.JCoups = AD_orcaS.JCoups
        AD_bobyqa.idx1Atoms = AD_orcaS.idx1Atoms
        AD_bobyqa.SParams = OrcaS_BOBYQA
        AD_bobyqa.method_save_files()
        print("\nThe data is saved to orcaS-BOBYQA.out file")
        print(f"{AD_bobyqa._file_orcaS}\n{AD_bobyqa.SParams}")
        print("  ========== End ==========")


if __name__ == "__main__":
    main()
