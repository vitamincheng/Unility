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
        default=15,
        help="Maximum of Combinations [default 16]",
    )

    parser.add_argument(
        "-l",
        "--limits",
        dest="limits",
        action="store",
        required=False,
        type=float,
        default=None,
        help="limits of border ppms [default 2 ppm in H, 20 ppm in C if .anmrrc is exists]",
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
    print("  ===== Loading data of peaks.npz from real spectra =====")
    peaks.method_read_file()
    peaks.method_print()

    # Average Directory of orcaS.out/orcaJ.out/orcaA.out and orcaS.BOBYQA
    AD_orcaS: AD_Normal = AD_Normal()
    if (AD_orcaS.Exist()):
        AD_orcaS.method_load_files()
        if isinstance(AD_orcaS.ChemicalShifts, dict):
            # convert dict to npt.NDArray to OrcaS
            OrcaS: npt.NDArray[np.float64] = np.array(
                list(AD_orcaS.ChemicalShifts.items()))
            print("\n  ===== Loading data OrcaS.out of Average Directory =====")
            print(f"{AD_orcaS._file_orcaS}\n{OrcaS}")
            in_SParams: list[int] = list(
                map(int, AD_orcaS.ChemicalShifts.keys()))
        else:
            print("  The OrcaS.out in Averaage Directory is not dict format !!!")
            print("  Exit and Close the program !!!")
            exit(0)

        from censo_ext.Tools.anmrfile import Anmr
        inAnmr: Anmr = Anmr()
        inAnmr.method_read_nucinfo()
        inAnmr.method_read_anmrrc()
        match(inAnmr.get_Anmrrc_Active()):
            case ["C"]:
                args.limits = 20
            case ["H"]:
                args.limits = 2

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
        Sim_ChemicalShifts: list[float] = []
        for x in unique_ChemEqvs_first_idx:
            index: tuple[npt.NDArray[np.intp], ...] = np.where(OrcaS.T[0] == x)
            Sim_ChemicalShifts.append(float(OrcaS.T[1][index][0]))

        sorted_ChemicalShifts: npt.NDArray[np.float64] = np.array(
            sorted(Sim_ChemicalShifts))
        real_Peaks: npt.NDArray[np.float64] = peaks.get_cIDs_center_peaks()
        print("\n  ===== Display the peaks of real spectra =====")
        print(f"Real Peaks : \n{real_Peaks.T}")
        print("\n  ===== Display the peaks of calculation spectra =====")
        print(f"Sorted ChemicalShifts : \n{sorted_ChemicalShifts}\n")

        # if len(sorted_ChemicalShifts) < len(real_Peaks.T):
        #    print(
        #        "  The numbers of real peaks are more than the numbers of simulation in normal")
        #    print("  Exit and Close the program !!!")
        #    exit(0)
        from censo_ext.Tools.utility import sub_numpy
        if len(sorted_ChemicalShifts) <= 4:
            print("  The numbers of sorted_data need more than four !!!")
            print("  Exit and Close the program !!!")
            exit(0)
        else:
            # nGroups: npt.NDArray[np.int64] = sub_numpy(
            #    sorted_ChemicalShifts, args.comb)
            if len(sorted_ChemicalShifts) >= args.comb:
                nGroups = sub_numpy(
                    sorted_ChemicalShifts, args.comb)
            else:
                nGroups = np.array([len(sorted_ChemicalShifts)])
        OrcaS_BOBYQA: npt.NDArray[np.float64] = np.insert(OrcaS, 2, 0, axis=1)

        print(f"  Sub-Groups : {nGroups}")

        counter: int = 1
        for idx, _ in enumerate(nGroups):

            start: np.int64 = np.sum(nGroups[:idx])+1
            if idx == 0:
                start = np.int64(0)
            end: np.int64 = np.sum(nGroups[:idx+1])

            Wait_Check_ChemicalShifts: npt.NDArray[np.float64] = sorted_ChemicalShifts[start:end+1]
            S_start: float = float(
                sorted_ChemicalShifts[start] - args.limits)
            S_end: float = float(sorted_ChemicalShifts[end] + args.limits)

            a: npt.NDArray[np.intp] = np.argwhere(real_Peaks[1] > S_start)
            b: npt.NDArray[np.intp] = np.argwhere(real_Peaks[1] < S_end)

            total_set: npt.NDArray[np.int64] = np.array(
                list(set(a.T[0]).intersection(set(b.T[0]))))

            if len(total_set) == 0:
                Wait_Check_Reals: npt.NDArray[np.float64] = np.array(
                    [])
            else:
                Wait_Check_Reals: npt.NDArray[np.float64] = np.array(
                    real_Peaks[1][total_set])

            # this process is for combination. it is more simple. And this is only for order numbers.
            if len(Wait_Check_Reals) >= len(Wait_Check_ChemicalShifts):
                Large = Wait_Check_Reals
                Small = Wait_Check_ChemicalShifts
            elif len(Wait_Check_ChemicalShifts) > len(Wait_Check_Reals):
                Large = Wait_Check_ChemicalShifts
                Small = Wait_Check_Reals
            else:
                raise ValueError("  BOBYQA_gen.py about 181 lines have error")
            # print(len(Large), len(Small))

            from itertools import combinations
            print("\n  ===== Combination of large nubmers of peaks =====")
            print(f"  The numbers of ({len(Large)}, {len(Small)})")

            combs: list[tuple[int, ...]] = list(combinations(
                [*range(0, len(Large))], len(Small)))  # type: ignore # nopep8
            result: list[float] = []
            print(f"  The numbers of Combinations : {len(combs)}")
            for comb in combs:
                result.append(cosine_similarity(
                    Large[list(comb)], Small))  # type: ignore # nopep8

            print(f"  The result of cosine_similarity : {np.array(result)}")

            print(f"\n  The best of the result : {max(result)}")
            args_combs: tuple[int, ...] = combs[result.index(max(result))]

            # print(len(Large), len(Small), args_combs)  # nopep8
            # ic(list(args_combs), Large, Small)
            if len(Wait_Check_Reals) >= len(Wait_Check_ChemicalShifts):
                # ordered_orcaS = Large[list(args_combs)]  # type: ignore # nopep8
                Wait_Check_Reals = Large[list(args_combs)]
                Wait_Check_ChemicalShifts = Small
            if len(Wait_Check_ChemicalShifts) > len(Wait_Check_Reals):
                # ordered_orcaS = Small[list(args_combs)]  # type: ignore # nopep8
                Wait_Check_ChemicalShifts = Large[list(args_combs)]
                Wait_Check_Reals = Small

                # Wait_Check_ChemicalShifts = Large[list(args_combs)]
            # exit(1001)
            # print(f"The best of the array : {ordered_orcaS}\n")
            print(f"  The best of the array (Real) : {Wait_Check_Reals}")
            print(
                f"  The best of the array (Cal.) : {Wait_Check_ChemicalShifts}\n")

            # ic(ordered_orcaS, OrcaS.T[1])
            for idx, x in enumerate(Wait_Check_ChemicalShifts):
                intp: tuple[npt.NDArray[np.intp], ...] = np.where(
                    OrcaS.T[1] == x)
                # intp: tuple[npt.NDArray[np.intp], ...] = np.where(
                #    OrcaS.T[1] == x)
                OrcaS_BOBYQA.T[1][intp] = Wait_Check_Reals[idx]

                if OrcaS_BOBYQA.T[2][intp] != 0:
                    print("repeated value ")
                    exit(0)
                else:
                    OrcaS_BOBYQA.T[2][intp] = counter
                counter += 1

            # for idx, x in enumerate(ordered_orcaS):
            #    intp: tuple[npt.NDArray[np.intp], ...] = np.where(
            #        real_Peaks[1] == x)
            #    # intp: tuple[npt.NDArray[np.intp], ...] = np.where(
            #    #    OrcaS.T[1] == x)
            #    OrcaS_BOBYQA.T[1][intp] = Wait_Check_Reals[idx]
            #    OrcaS_BOBYQA.T[2][intp] = counter
            #    counter += 1

        # ic(OrcaS_BOBYQA)
        AD_bobyqa: AD_BOBYQA = AD_BOBYQA()
        AD_bobyqa.JCoups = AD_orcaS.JCoups
        AD_bobyqa.idx1Atoms = AD_orcaS.idx1Atoms
        AD_bobyqa.ChemicalShifts = OrcaS_BOBYQA
        AD_bobyqa.method_save_files()
        print("\nThe data is saved to orcaS-BOBYQA.out file")
        print(f"{AD_bobyqa._file_orcaS}\n{AD_bobyqa.ChemicalShifts}")
        print("  ========== End ==========")


if __name__ == "__main__":
    main()
