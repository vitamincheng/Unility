#!/usr/bin/env python
from censo_ext.Tools.anmrfile import AD_BOBYQA, AD_Normal
from censo_ext.Tools.utility import AtomID, print_arguments
from icecream import ic
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
        default=100,
        help="Maximum of Combinations [default 100]",
    )

    parser.add_argument(
        "-start",
        "--start",
        dest="start",
        action="store",
        required=False,
        type=float,
        help="start ppm of Chemical shift",
    )

    parser.add_argument(
        "-end",
        "--end",
        dest="end",
        action="store",
        required=False,
        type=float,
        help="end ppm of Chemical shift",
    )

    parser.add_argument(
        "-d",
        "--del",
        dest="delete",
        action="store",
        required=False,
        type=int,
        nargs="+",
        default=None,
        help="Under Calculation, the number of neglect atoms in orcaS.out",
    )

    parser.add_argument(
        "--index",
        dest="index",
        action="store",
        required=False,
        type=int,
        nargs="+",
        default=None,
        help="index of the peaks in peaks.npz [default None]",
    )

    return parser.parse_args()


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
            print(f"{AD_orcaS._orcaS}\n{OrcaS}")
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
        Anmrrc_Active: list[str] = inAnmr.get_Anmrrc_Active()

        if args.start and args.end and not args.index:
            pass
        else:
            if Anmrrc_Active == ['H']:
                limits = 0.5
            elif Anmrrc_Active == ['C']:
                limits = 10
            else:
                print(" Active element of anmrrc is not H or C ")
                exit(0)

        ChemEqvs: dict[AtomID, list[AtomID]] = {key: value for key, value in
                                                inAnmr.NeighborChemEqvs.items()
                                                if key in in_SParams}
        Sorted_ChemEqvs: list[list[AtomID]] = list(
            sorted(value) for value in ChemEqvs.values())

        unique_ChemEqvs: list[list[AtomID]] = []

        for Sorted_ChemEqv in Sorted_ChemEqvs:
            if Sorted_ChemEqv not in unique_ChemEqvs:
                unique_ChemEqvs.append(Sorted_ChemEqv)

        unique_ChemEqvs_first_idx: list[AtomID] = []
        for item in unique_ChemEqvs:
            unique_ChemEqvs_first_idx.append(item[0])

        # Normal Sim_SParams/sorted_SParams is more than real known peaks
        Sim_CS: list[float] = []
        for x in unique_ChemEqvs_first_idx:
            _idx: tuple[npt.NDArray[np.intp], ...] = np.where(
                OrcaS.T[0] == x)
            Sim_CS.append(float(OrcaS.T[1][_idx][0]))

        sorted_CS: npt.NDArray[np.float64] = np.array(
            sorted(Sim_CS))
        real_Peaks: npt.NDArray[np.float64] = peaks.get_cIDs_center_peaks()
        print("\n  ===== Display the peaks of real spectra =====")
        print(f"The numbers of Real Peaks : {len(real_Peaks.T)}")
        print(f"Real Peaks : \n{real_Peaks.T}")
        print("\n  ===== Display the peaks of calculation spectra =====")
        print(f"The numbers of Sorted ChemicalShifts : {len(sorted_CS)}")
        print(f"Sorted ChemicalShifts : \n{sorted_CS}\n")

        OrcaS_BOBYQA: npt.NDArray[np.float64] = np.insert(
            OrcaS, 2, 0, axis=1)

        # Wait_Check_CS: npt.NDArray[np.float64] = sorted_CS
        Wait_Check_Reals: npt.NDArray[np.float64] = real_Peaks[1]

        if not args.index:
            if args.start and args.end:
                start = args.start
                end = args.end
            else:
                start = Wait_Check_Reals.min() - limits  # type: ignore
                end = Wait_Check_Reals.max() + limits  # type: ignore

            args_start = sorted_CS > start
            args_end = sorted_CS < end
            args_intersection = np.logical_and(args_start, args_end)
            Wait_Check_CS = sorted_CS[args_intersection]

            if args.delete is not None:
                print("  Activated Delete Atoms : ")
                for x in args.delete:
                    y = OrcaS[np.where(OrcaS.T[0] == x)[0]][0][1]
                    print(f"  {x}  {y}")
                    z = np.where(Wait_Check_CS == y)[0]
                    Wait_Check_CS = np.delete(Wait_Check_CS, z)
            print(
                f"The numbers of be Checked ChemicalShifts : {len(Wait_Check_CS)}")
            print(f"be Checked ChemicalShifts : \n{Wait_Check_CS}\n")

        else:
            list_x: list[int] = []
            for x in args.index:
                a = np.argwhere(OrcaS.T[0] == x)[0][0]
                list_x.append(int(a))
            Wait_Check_CS: npt.NDArray = OrcaS[list_x].T[1]
            Wait_Check_CS.sort()

        # this process is for combination. it is more simple. And this is only for order numbers.
        # Normally the distance of Wait_Check_Reals is more width than Wait_Check_CS, so peaks is more
        # If the peaks of Wait_Check_Reals is not including in. it will less.
        if len(Wait_Check_Reals) >= len(Wait_Check_CS):
            Large = Wait_Check_Reals
            Small = Wait_Check_CS
        elif len(Wait_Check_CS) > len(Wait_Check_Reals):
            Large = Wait_Check_CS
            Small = Wait_Check_Reals
        else:
            ic()
            raise ValueError("  BOBYQA_gen.py have error")

        from itertools import combinations
        print("\n  ===== Combination of large nubmers of peaks =====")
        print(f"  The numbers of ({len(Large)}, {len(Small)})")

        combs: list[tuple[int, ...]] = list(combinations(
            [*range(0, len(Large))], len(Small)))  # type: ignore # nopep8
        result: list[float] = []
        print(f"  The numbers of Combinations : {len(combs)}")
        for comb in combs:
            distance = Large[list(comb)]-Small
            result.append(np.sum(np.square(distance)))

        print(f"  The result of Euclidean Distance : {np.array(result)}")

        print(f"\n  The best of the result : {min(result)}")
        args_combs: tuple[int, ...] = combs[result.index(min(result))]

        if len(Wait_Check_Reals) >= len(Wait_Check_CS):
            Wait_Check_Reals = Large[list(args_combs)]
            Wait_Check_CS = Small
        if len(Wait_Check_CS) > len(Wait_Check_Reals):
            Wait_Check_CS = Large[list(args_combs)]
            Wait_Check_Reals = Small

        print(f"  The best of the array (Real) : {Wait_Check_Reals}")
        print(
            f"  The best of the array (Cal.) : {Wait_Check_CS}\n")

        counter = 1
        for idx, x in enumerate(Wait_Check_CS):
            intp: tuple[npt.NDArray[np.intp], ...] = np.where(
                OrcaS.T[1] == x)
            OrcaS_BOBYQA.T[1][intp] = Wait_Check_Reals[idx]

            if OrcaS_BOBYQA.T[2][intp].all() != 0:
                print("repeated value ")
                exit(0)
            else:
                OrcaS_BOBYQA.T[2][intp] = counter
            counter += 1

        AD_bobyqa: AD_BOBYQA = AD_BOBYQA()
        AD_bobyqa.JCoups = AD_orcaS.JCoups
        AD_bobyqa.Element = AD_orcaS.Element
        AD_bobyqa.ChemicalShifts = OrcaS_BOBYQA
        AD_bobyqa.method_save_files()
        print("\nThe data is saved to orcaS-BOBYQA.out file")
        print(f"{AD_bobyqa._orcaS}\n{AD_bobyqa.ChemicalShifts}")
        print("  ========== End ==========")


if __name__ == "__main__":
    main()
