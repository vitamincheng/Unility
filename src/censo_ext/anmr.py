#!/usr/bin/env python
import copy
import numpy as np
import numpy.typing as npt
import argparse
from icecream import ic
from pathlib import Path

descr = """
________________________________________________________________________________
| anmr.py  
| Usages   : anmr.py [options]
| [options]
| Output   : -o output file [default output.dat]
| Dir      : -D the directory of input files(CONF) [default .]
| mf       : -mf magnetic frequency of scan nmr [default 500.0]
| lw       : -lw line width of scan nmr [1.0 for H, 20 for C]
| auto     : -auto --auto automated to adjust the threshold of J and AB quartet [default False]  
| average  : -av load the average folder data to plot spectra [default False] 
| thr      : -t -thr threshold of coupling constant (J) [default 0.30]
| thrtab   : -tab -thrab threshold of AB quartet (J/d chemical shift) [default 0.025]
| tbpent   : -tb threshold of AB quartet bond pententration distance [default 4]
| mss      : -mss max of spin numbers [default 10]
|            if your computer have slowly CPU, try to use mss 4 or 5. 
| Cutoff   : --cutoff cutoff limits in quantum chemistry in nmr [default 0.001]
| Start    : -start start ppm of plotting spectra
| End      : -end end ppm of plotting spectra
| ref      : reference standard - see .anmrrc file
| BOBYQA   : -b --bobyqa BOBYQA mode (only use numpy) [default False]
| JSON     : -j --json Read the raw data of every single peak [if is -1(All)]
| ascal    : -ascal chemical shift scaling a if the reference is absent [pending]
| bscal    : -bcsal chemical shift scaling b if the reference is absent [pending]
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
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.dat",
        help="Provide output_file name [default output.dat]",
    )

    parser.add_argument(
        "-D",
        "--dir",
        dest="dir",
        action="store",
        required=False,
        default=".",
        help="Provide the directory name [default .]",
    )

    parser.add_argument(
        "-mf",
        "--magnfreq",
        dest="mf",
        action="store",
        type=float,
        required=False,
        default=500.0,
        help="magnetic frequency of scan nmr [default 500.0]",
    )
    parser.add_argument(
        "-lw",
        "-linewidth",
        dest="lw",
        action="store",
        type=float,
        required=False,
        help="line width of scan nmr [default 1.0 for H, 20.0 for C]",
    )

    parser.add_argument(
        "-ascal",
        dest="ascal",
        action="store",
        type=float,
        required=False,
        help="chemical shift scaling a, "
        "For   "
        "For   ",
    )

    parser.add_argument(
        "-bscal",
        dest="bscal",
        action="store",
        type=float,
        required=False,
        help="chemical shift scaling b, "
        "For  "
        "For  ",
    )

    parser.add_argument(
        "-t",
        "--thr",
        dest="thr",
        action="store",
        type=float,
        required=False,
        help="threshold of J coupling constant Hz [default 0.30]",
    )

    parser.add_argument(
        "-tab",
        "--thrab",
        dest="thrab",
        action="store",
        type=float,
        required=False,
        default=0.025,
        help="threshold of AB quartet (J/chemical shift) [default 0.025]",
    )

    parser.add_argument(
        "-tb",
        "--tbpent",
        dest="tb",
        action="store",
        type=int,
        required=False,
        default=4,
        help="threshold of AB quartet bond pententration distance [default 4]",
    )

    parser.add_argument(
        "-cut",
        "--cutoff",
        dest="cutoff",
        action="store",
        type=float,
        default=0.001,
        help="cutoff limits in quantum chemistry in nmr [default 0.001]",
    )

    parser.add_argument(
        "-start",
        "--start",
        dest="start",
        action="store",
        type=float,
        required=False,
        help="start ppm of plotting spectra",
    )

    parser.add_argument(
        "-end",
        "--end",
        dest="end",
        action="store",
        type=float,
        required=False,
        help="end ppm of plotting spectra",
    )

    parser.add_argument(
        "-av",
        "--average",
        dest="average",
        action="store_true",
        help="Load the averager foler data to plot spectra when use auto argument [default False]",
    )

    parser.add_argument(
        "-b",
        "--bobyqa",
        dest="bobyqa",
        action="store_true",
        help="BOBYQA mode (only use numpy) [default False]",
    )

    parser.add_argument(
        "-j",
        "--json",
        dest="json",
        action="store",
        type=int,
        nargs="+",
        required=False,
        help="Read the jason file of raw single peak [if is -1(All)]",
    )

    parser.add_argument(
        "-mss",
        "--mss",
        dest="mss",
        action="store",
        type=int,
        required=False,
        default=10,
        help="max of spin numbers [default 10]",
    )

    parser.add_argument(
        "-auto",
        "--auto",
        dest="auto",
        action="store_true",
        help="Automated to adjust the threshold of AB quartet [defalut False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> npt.NDArray[np.float64]:

    if args == argparse.Namespace():
        args = cml("")

    Directory: Path = Path("")
    if args.dir:
        Directory = Path(args.dir)

    # Import necessary modules for processing
    from censo_ext.Tools.anmrfile import Anmr

    # Create an Anmr object and read required files
    inAnmr: Anmr = Anmr(Directory)
    inAnmr.method_read_anmrrc()
    inAnmr.method_read_nucinfo()

    # Handle average data loading if specified
    if inAnmr.get_avg_orcaSJ_Exist() and args.average:
        inAnmr.method_load_avg_orcaSJ(args=args)
    else:
        # Process all ORCA files and generate average data
        inAnmr.method_read_enso()
        inAnmr.method_read_folder_orcaSJ()
        inAnmr.method_update_equiv_orcaSJ()
        inAnmr.method_avg_orcaSJ()
        inAnmr.method_save_avg_orcaSJ()

    # Extract spin parameters and coupling constants
    inSParams: npt.NDArray[np.float64] = np.array(
        list(inAnmr.avg_orcaSJ.SParams.values()))*args.mf

    # Get coupling constants and atom information
    inJCoups: npt.NDArray[np.float64] = np.array(inAnmr.avg_orcaSJ.JCoups)
    in_idx1Atoms: dict[int, str] = inAnmr.avg_orcaSJ.idx1Atoms

    # Initialize variables for processing based on active nuclear element
    dpi: int | None = None
    Active_range: int | None = None
    inHydrogen: list[int] = []
    inFile: Path = Path("crest_conformers.xyz")

    # Process different nuclear elements (C or H)
    for idx, Active in enumerate(inAnmr.get_Anmr_Active()):
        if idx == 0:
            if Active == 'C':
                # Carbon processing - read molecular structure
                from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order
                mol, neighbors, bond_order = read_mol_neighbors_bond_order(
                    inAnmr.get_Directory() / inFile)

                # For C/CH/CH2/CH3 from 0 1 2 3 to 1 2 3 4 for Carbon spectra
                inHydrogen = [value+1 for value in bond_order.values()]

                Active_range = 200
                dpi = 500
                if not args.lw:
                    args.lw = 20
                if not args.thr:
                    args.thr = args.lw * 0.3
                inJCoups = np.zeros_like(inJCoups)
            elif Active == 'H':
                # Hydrogen processing - identify equivalent hydrogens
                idx1_nMagEqvHydrogens: dict[int, int] = {}
                for key in inAnmr.avg_orcaSJ.idx1Atoms.keys():
                    idx1_nMagEqvHydrogens[key] = inAnmr.nMagnetEqvs[key]
                for y in inAnmr.get_idx1_acid_atoms_NoShow_RemoveH(
                        inAnmr.get_Directory()/inFile):
                    if y in idx1_nMagEqvHydrogens.keys():
                        del idx1_nMagEqvHydrogens[y]
                        del in_idx1Atoms[y]

                inHydrogen = [
                    value for value in idx1_nMagEqvHydrogens.values()]

                Active_range = 10
                dpi = 10000
                if not args.lw:
                    args.lw = 1
                if not args.thr:
                    args.thr = args.lw * 0.3
            else:
                print("  Other Active Nuclear element, waiting to build")
                print("  Exit and Close the program !!!")
                ic()
                exit(0)
        elif idx >= 1:
            print("  Only for ONE Active Nuclear element, waiting to build")
            print("  Exit and Close the program !!!")
            ic()
            exit(0)

    # Initialize variables for AB quartet detection and processing
    idx0_ab_group_sets: list[set[int]] = []
    mat_filter_multi: npt.NDArray[np.int64] = np.array([])

    # Main processing loop for identifying and categorizing spin systems
    if not args.json and len(inAnmr.get_Anmr_Active()) == 1 and inAnmr.get_Anmr_Active()[0] == 'H':
        inJCoups_origin: npt.NDArray[np.float64] = copy.deepcopy(inJCoups)

        while (1):

            # the numbers of mss is low, it will be computated as AB quartet
            if len(inSParams*inHydrogen) <= args.mss:

                inJCoups = copy.deepcopy(inJCoups_origin)

                # Step 1: Filter out small coupling constants based on threshold
                # Delete Too Small JCoups J = args.lw*(-0.3) ~ args.lw*(0.3) use matrix Filter
                # 1: keep and 0: neglect
                # reset all inJCoups
                mat_filter_low_factor: npt.NDArray[np.int64] = np.zeros(
                    (inSParams.size, inSParams.size), dtype=np.int64)
                mat_filter_low_factor = (
                    np.abs(inJCoups) > args.thr).astype(np.int64)
                inJCoups[np.logical_not(mat_filter_low_factor)] = 0
                del mat_filter_low_factor

                # Step 2: add more one hydrogen in the same position as the first hydrogen
                idx_repeat = (np.array(inHydrogen)-1).nonzero()[0]
                for idx in idx_repeat[::-1]:
                    for idy in range(1, inHydrogen[idx], 1):
                        inSParams = np.insert(inSParams, idx, inSParams[idx])
                        inJCoups = np.insert(inJCoups, idx, inJCoups[idx], axis=1)  # nopep8
                        inJCoups = np.insert(inJCoups, idx, inJCoups[idx], axis=0)  # nopep8
                inHydrogen = [1]*len(inSParams)
                mat_filter_low_factor = np.ones(
                    (inSParams.size, inSParams.size), dtype=np.int64)
                np.fill_diagonal(mat_filter_low_factor, 0)
                mat_filter_ab_quartet: npt.NDArray[np.int64] = copy.deepcopy(
                    mat_filter_low_factor)
            else:
                # Step 1: Filter out small coupling constants based on threshold
                # Delete Too Small JCoups J = args.lw*(-0.3) ~ args.lw*(0.3) use matrix Filter
                # 1: keep and 0: neglect
                mat_filter_low_factor: npt.NDArray[np.int64] = np.zeros(
                    (inSParams.size, inSParams.size), dtype=np.int64)
                inJCoups = copy.deepcopy(inJCoups_origin)
                mat_filter_low_factor = (
                    np.abs(inJCoups) > args.thr).astype(np.int64)
                inJCoups[np.logical_not(mat_filter_low_factor)] = 0
                # mat_filter_low_factor = (np.abs(inJCoups_origin) > args.thr)*1
                # inJCoups = copy.deepcopy(inJCoups_origin)
                # inJCoups *= mat_filter_low_factor

                # Step 2: Identify potential AB quartet systems
                # np.fill_diagonal(mat_filter_ab_quartet, 0)
                mat_filter_ab_quartet: npt.NDArray[np.int64] = np.zeros(
                    (inSParams.size, inSParams.size), dtype=np.int64)
                import math
                for idx, x in enumerate(inSParams):
                    for idy, y in enumerate(inSParams):
                        if idx == idy:
                            mat_filter_ab_quartet[idx][idy] = 0
                        else:
                            # Check if two chemical shifts are very close (AB quartet condition)
                            # If the JCoups is negative, the pwaks will prioritize to normal QM calculation
                            # if two chemical shift is very close, will perform AB quartet
                            # if x-y == 0 the Ratio_J_Hz will crash
                            # if the JCoups is negative, will prioritize to use normal QM calculation.
                            # if the peaks which is not AB Quartet and have positive nubmers will use Multiplet (nmrsim)
                            if (math.fabs(x-y) < 0.0005 and mat_filter_low_factor[idx][idy] == 1) or (inJCoups[idx][idy] <= -args.thr):
                                mat_filter_ab_quartet[idx][idy] = 1
                            else:
                                if math.fabs(x-y) < 0.0005:
                                    Ratio_J_Hz = 10000
                                else:
                                    Ratio_J_Hz = math.fabs(
                                        inJCoups[idx][idy]/(math.fabs(x-y)))
                                if Ratio_J_Hz < args.thrab:
                                    mat_filter_ab_quartet[idx][idy] = 0
                                elif Ratio_J_Hz >= args.thrab and mat_filter_low_factor[idx][idy] == 1:
                                    mat_filter_ab_quartet[idx][idy] = 1
                                else:
                                    if mat_filter_low_factor[idx][idy] == 1:
                                        ic(idx, x, idy, y,
                                           "  Something wrong in your SParams !!!")
                                        ic()
                                        raise ValueError(
                                            f"{idx} {x} {idy} {y} was not found or is a directory")

            # ic(mat_filter_ab_quartet)
            # Calculate which couplings are NOT part of AB quartets (multiplets)
            mat_filter_multi = mat_filter_low_factor - mat_filter_ab_quartet
            idx0_ab_connect: list[list[int | set[int]]] = []
            for idx0, x in enumerate(mat_filter_ab_quartet):
                group: set[int] = set(np.array(x*(idx0+1)).nonzero()[0].tolist())  # return arg # nopep8 #idx0+1 is only for nozero
                group.add(idx0)
                idx0_ab_connect.append([idx0, group])

            # Merge overlapping spin system groups
            idx0_ab_group_sets = []
            for _, x in idx0_ab_connect:
                if len(x) == 1:                     # type: ignore
                    idx0_ab_group_sets.append(x)     # type: ignore
                elif len(x) > 1:                    # type: ignore
                    group: set[int] = x             # type: ignore
                    loop: bool = True
                    bond_penetration: int = 1
                    while loop:
                        loop = False
                        for idy, y in enumerate(group):
                            if (not group.issuperset(idx0_ab_connect[y][1])) and bond_penetration <= args.tb:  # type: ignore # nopep8
                                loop = True
                                group = group.union(idx0_ab_connect[y][1])  # type: ignore # nopep8
                            bond_penetration += 1
                    idx0_ab_group_sets.append(group)
                else:
                    print("  Something wrong in ab_group_set")
                    ic()
                    raise ValueError("  idx0_ab_group_sets is bugs !!!")

            if len(inSParams*inHydrogen) <= args.mss:
                pass
            else:
                # Handle CH3 equivalent groups manually (symmetry considerations)
                # So if chemical shift in AB quartet region need to move to multiplet
                list_Equivalent3: list[int] = []
                for key in inAnmr.nMagnetEqvs.keys():
                    if inAnmr.nMagnetEqvs[key] == 3:
                        for idy, y in enumerate(in_idx1Atoms):
                            if y == min(inAnmr.NeighborMangetEqvs[key]):
                                list_Equivalent3.append(idy)
                set_Equivalent3: set[int] = set(list_Equivalent3)

                # Adjust groups to account for equivalent protons
                # Equivalent3 is idx0 numbers
                for idx, group_set in enumerate(idx0_ab_group_sets):
                    set_move: set[int] = group_set.intersection(
                        set_Equivalent3)
                    if not len(set_move) == 0:
                        idx0_ab_group_sets[idx] = set(group_set).difference(
                            set_move).union(set([idx]))
                        set_move = set_move.difference(set([idx]))
                    for idy, y in enumerate(set_move):
                        mat_filter_multi[idx][y] = 1
            # ic(idx0_ab_group_sets)
            print(" ===== Processing =====")
            print(f" threshold of JCoupling  : {args.thr:>3.5f}")
            print(f" threshold of AB quartet : {args.thrab:>3.5f}")
            print("  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")

            max_len_AB: int = 0
            for idx, group_set in enumerate(idx0_ab_group_sets):
                mat_multi_x_idx: list[int] = [
                    idx0_set*x for x, idx0_set in enumerate(mat_filter_multi[idx].tolist())if idx0_set != 0]
                print(f'{(idx+1):>5d}{len(group_set):>5d}', {a+1 for a in group_set}, set(
                    a+1 for a in mat_multi_x_idx).difference({a+1 for a in group_set}))
                if len(group_set) > max_len_AB:
                    max_len_AB = len(group_set)

            # Adjust thresholds automatically if maximum spin system exceeds limit
            if (max_len_AB > args.mss):
                if args.auto:
                    args.thrab = args.thrab + 0.0025
                    args.thr = args.thr + 0.030
                else:
                    print("  Need to tidy the nspins of AB quartet and use cmd -thrab ")
                    print("  Exit and Close the program !!!")
                    exit(0)
            else:
                break

        # Additional processing for AB quartets with identical chemical shifts
        # AB quartet if more than two peaks of AB quartet, added closed peaks (not in AB quartet)
        print(" ===== Modification AB quartet =====")
        for idx0, idx0_ab_group_set in enumerate(idx0_ab_group_sets):
            idx0_ab_group: list[int] = list(idx0_ab_group_set)
            if len(set(list(inSParams[idx0_ab_group]))) == 1:
                from censo_ext.Tools.spectra import find_nearest
                _, Move_idx0 = find_nearest(list(inSParams[idx0_ab_group]),
                                            inSParams[idx0_ab_group[0]])
                a = np.argwhere(
                    inSParams[:] == inSParams[int(Move_idx0)])
                idx0_ab_group_sets[idx0] = idx0_ab_group_set.union(
                    set(int(idx) for idx in a[0]))

        # Display the parameter of Full Spectra
        for idx0, idx0_ab_group_set in enumerate(idx0_ab_group_sets):
            idx0_ab_group: list[int] = list(idx0_ab_group_set)
            mat_multi_x_idx0: list[int] = [
                idx0_set*x for x, idx0_set in enumerate(mat_filter_multi[idx0].tolist())if idx0_set != 0]
            print(f'{(idx0+1):>5d}{len(idx0_ab_group):>5d}', {a+1 for a in idx0_ab_group}, set(
                a+1 for a in mat_multi_x_idx0).difference({a+1 for a in idx0_ab_group}))
        print(" Use this parameter to calculate the Full Spectra")

    # Low level QM model
    # see https://nmrsim.readthedocs.io/en/latest/index.html
    #
    import censo_ext.Tools.qm as qm
    idx0_peaks_range: list[int] = []
    Res_peaks: list[list[tuple[float, float]]] = []
    if not args.json and len(inAnmr.get_Anmr_Active()) == 1 and inAnmr.get_Anmr_Active()[0] == 'H':
        print("")
        print(" ===== Processing =====")
        print(" the group of calculate spectra :", len(idx0_ab_group_sets))
        print("  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")
        # Res_peaks: list[list[tuple[float, float]]] = []
        if len(inSParams*inHydrogen) <= args.mss:
            idx0_ab_group = list(idx0_ab_group_sets[0])
            v: npt.NDArray[np.float64] = inSParams[idx0_ab_group]
            J: npt.NDArray[np.float64] = inJCoups[idx0_ab_group].T[idx0_ab_group]
            QM_Base: list[tuple[float, float]] = qm.qm_full(v=list(
                v), J=J, nIntergals=len(inHydrogen), args=args)
            from nmrsim.math import normalize_peaklist
            QM_Multiplet = normalize_peaklist(QM_Base, len(inHydrogen))
            Res_peaks.append(QM_Multiplet)

        else:
            for idx0, idx0_ab_group_set in enumerate(idx0_ab_group_sets):
                mat_multi_idx0: list[int] = mat_filter_multi[idx0].astype(
                    int).tolist()
                idx0_ab_group: list[int] = list(idx0_ab_group_set)
                idx1_ab_group: set[int] = {a+1 for a in idx0_ab_group_set}
                mat_multi_x_idx0: list[int] = [
                    idx0_set*x for x, idx0_set in enumerate(mat_multi_idx0)if idx0_set != 0]
                print(f'{(idx0+1):>5d}{len(idx0_ab_group):>5d}', f'{idx1_ab_group}', set(
                    a+1 for a in mat_multi_x_idx0).difference(idx1_ab_group))

                v: npt.NDArray[np.float64] = inSParams[idx0_ab_group]
                J: npt.NDArray[np.float64] = inJCoups[idx0_ab_group].T[idx0_ab_group]

                QM_Base: list[tuple[float, float]] = qm.qm_base(v=list(
                    v), J=J, nIntergals=inHydrogen[idx0_ab_group.index(idx0)], idx0_nspins=idx0_ab_group.index(idx0), args=args)

                QM_Multiplet: list[tuple[float, float]] = []
                for z in QM_Base:
                    multiplicity: list[int] = list(
                        set(mat_multi_x_idx0).difference(idx0_ab_group_set))
                    inJ: list[tuple[float, int]] = []
                    for a in multiplicity:
                        if np.fabs(inSParams[idx0]-inSParams[a]) > 0.1:
                            inJ.append((inJCoups[idx0][a], inHydrogen[a]))

                    if len(inJ) >= 1:
                        tmp: npt.NDArray[np.float64] = np.array(
                            qm.qm_multiplet(z[0], nIntergals=1, J=inJ))
                        tmp.T[1] *= z[1]
                        QM_Multiplet += tmp.tolist()
                    elif len(inJ) == 0:
                        QM_Multiplet = QM_Base
                    else:
                        print("Something wrong")
                        print("  Exit and Close the program !!!")
                        ic()
                        exit(1)
                from nmrsim.math import normalize_peaklist
                QM_Multiplet = normalize_peaklist(
                    QM_Multiplet, inHydrogen[idx0])

                if len(Res_peaks) == 0 and len(QM_Multiplet) == 0:
                    pass
                else:
                    Res_peaks.append(QM_Multiplet)

        import json
        # import pickle
        with open(str(inAnmr.get_Directory()/Path("peaks.json")), "w") as final:
            json.dump(Res_peaks, final)

        idx0_peaks_range = [*range(len(Res_peaks))]

    elif not args.json and len(inAnmr.get_Anmr_Active()) == 1 and inAnmr.get_Anmr_Active()[0] == 'C':
        for idx, ppm in enumerate(inSParams):
            dat: list = []
            dat.append((float(ppm*(-1)), float(inHydrogen[idx])))
            Res_peaks.append(dat)

        import json
        # import pickle
        with open(str(inAnmr.get_Directory()/Path("peaks.json")), "w") as final:
            json.dump(Res_peaks, final)
    else:
        import json
        # import pickle
        with open(str(inAnmr.get_Directory()/Path("peaks.json")), "r") as json_file:
            Res_peaks = json.load(json_file)

        if args.json[0] == -1:
            idx0_peaks_range = [*range(len(Res_peaks))]
        else:
            idx0_peaks_range = args.json

    Final_peaks: list[tuple[float, float]] = []

    if inAnmr.get_Anmr_Active()[0] == 'H':
        for idx, peak in enumerate(Res_peaks):
            if idx in idx0_peaks_range:
                Final_peaks += peak
    elif inAnmr.get_Anmr_Active()[0] == 'C':
        for idx, peak in enumerate(Res_peaks):
            Final_peaks += peak
    else:
        print("  Something Wrong in your get_anmr_Active()")
        ic()
        exit(0)

    print("")
    print(" ===== Processing to plotting spectra =====")
    print(" Wait a minutes ...")
    args.out = str(inAnmr.get_Directory()/Path(args.out))

    if dpi and Active_range:
        np_dat: npt.NDArray[np.float64] = qm.print_plot(
            Final_peaks, dpi, nIntergals=1, args=args, Active_range=Active_range)
    else:
        print("  dpi and Active_range is wrong")
        print("  Exit and Close the program !!!")
        ic()
        exit(0)

    if not args.bobyqa:
        print(f" the spectra is saved to : {args.out}")

    print(" All done ...")
    return np_dat


if __name__ == "__main__":
    main()
