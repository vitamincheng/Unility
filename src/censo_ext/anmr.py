#!/usr/bin/env python
import os
import copy
from typing import Union
import numpy as np
import numpy.typing as npt
import argparse
from icecream import ic
from pathlib import Path

descr = """
________________________________________________________________________________
|                                          [10.16.2024] vitamin.cheng@gmail.com
| anmr.py  
| Usages  : anmr.py [options]
| [options]
| Output   : -o output file [default output.dat]
| Dir      : -d the directory of input files(CONF) [default .]
| mf       : -mf magnetic frequency of scan nmr [default 500.0]
| lw       : -lw line width of scan nmr [1.0 for H, 20 for C]
| auto     : -auto --auto automated to adjust the threshold of J and AB quartet   
| average  : -av load the average folder data to plot spectra  
| thr      : -t -thr threshold of coupling constant (J) [default 0.30]
| thrtab   : -tab -thrab threshold of AB quartet (J/d chemical shift) [default 0.025]
| tbpent   : -tb threshold of AB quartet bond pententration distance [default 4]
| mss      : -mss max of spin numbers [default 9]
|            if your computer have slowly CPU, try to use mss 4 or 5. 
| Cutoff   : --cutoff cutoff limits in quantum chemistry in nmr [default 0.001]
| Start    : -start start ppm of plotting spectra
| End      : -end end ppm of plotting spectra
| ref      : reference standard - see .anmrrc file
| BOBYQA   : -b --bobyqa BOBYQA mode (only use numpy) [default False]
| JSON     : -j --json Read the raw data of every single peak 
| ascal    : -ascal chemical shift scaling a if the reference is absent [pending]
| bscal    : -bcsal chemical shift scaling b if the reference is absent [pending]
| Package  : Tools / nmrsim Library
| Module   : anmrfile.py / ml4nmr.py / qm.py
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
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.dat",
        help="Provide output_file name [default output.dat]",
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
        default=False,
        help="Load the averager foler data to plot spectra when use auto argument",
    )

    parser.add_argument(
        "-show",
        "--show",
        dest="show",
        action="store_true",
        default=False,
        help="Show the spectra",
    )

    parser.add_argument(
        "-b",
        "--bobyqa",
        dest="bobyqa",
        action="store_true",
        help="BOBYQA mode (only use numpy)",
    )

    parser.add_argument(
        "-j",
        "--json",
        dest="json",
        action="store",
        type=int,
        nargs="+",
        required=False,
        help="Read the jason file of raw single peak",
    )

    parser.add_argument(
        "-mss",
        "--mss",
        dest="mss",
        action="store",
        type=int,
        required=False,
        default=9,
        help="max of spin numbers [default 9]",
    )

    parser.add_argument(
        "-auto",
        "--auto",
        dest="auto",
        action="store_true",
        default=False,
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

    from censo_ext.Tools.anmrfile import Anmr
    inAnmr: Anmr = Anmr(Directory)
    inAnmr.method_read_anmrrc()
    inAnmr.method_read_nucinfo()

    if inAnmr.get_avg_orcaSJ_Exist() and args.average:
        inAnmr.method_load_avg_orcaSJ(args=args)
    else:
        inAnmr.method_read_enso()
        inAnmr.method_read_folder_orcaSJ()
        inAnmr.method_update_equiv_orcaSJ()
        inAnmr.method_avg_orcaSJ()
        inAnmr.method_save_avg_orcaSJ()

    inSParams: npt.NDArray[np.float64] = np.array(
        list(inAnmr.avg_orcaSJ.SParams.values()))*args.mf

    # idxinSParams:npt.NDArray[np.float64]  = np.array(list(inAnmr.Average_orcaSJ.orcaSParams.keys()))
    # np.savetxt("inSParams.out",np.array(np.stack([idxinSParams,inSParams]).T), fmt=' %4.0f   %10.6f')

    inJCoups: npt.NDArray[np.float64] = np.array(inAnmr.avg_orcaSJ.JCoups)
    in_idxAtoms: dict[int, str] = inAnmr.avg_orcaSJ.idxAtoms
    dpi: int | None = None
    Active_range: int | None = None
    inHydrogen: list[int] = []
    in_File: Path = Path("crest_conformers.xyz")
    # ic(inAnmr.Directory)
    for idx, x in enumerate(inAnmr.get_Anmr_Active()):
        if idx == 0:
            if x == 'C':
                from censo_ext.Tools.ml4nmr import read_mol_neighbors_bond_order
                mol, neighbors, bond_order = read_mol_neighbors_bond_order(
                    inAnmr.get_Directory() / in_File)
                inHydrogen = [value+1 for value in bond_order.values()]
                Active_range = 200
                dpi = 500
                if args.lw == None:
                    args.lw = 20
                if not args.thr:
                    args.thr = args.lw * 0.3
                inJCoups = np.zeros_like(inJCoups)
            elif x == 'H':
                inHydrogenDict: dict[int, int] = {}
                for key in inAnmr.avg_orcaSJ.idxAtoms.keys():
                    inHydrogenDict[key] = len(inAnmr.nucinfo[1][key-1][2])
                for x in inAnmr.get_idx1_acid_atoms_NoShow_RemoveH(
                        inAnmr.get_Directory()/in_File):
                    if x in inHydrogenDict.keys():
                        del inHydrogenDict[x]
                        del in_idxAtoms[x]
                inHydrogen = [value for value in inHydrogenDict.values()]
                Active_range = 10
                dpi = 10000
                if args.lw == None:
                    args.lw = 1
                if args.thr == None:
                    args.thr = args.lw * 0.3
            else:
                print("Other Active Nuclear element, waiting to build")
                ic()
                exit(0)
        elif idx >= 1:
            print("only for ONE Active Nuclear element, waiting to build")
            ic()
            exit(0)
    idx0_ab_group_set: list[set[int]] = []
    mat_filter_multi: npt.NDArray[np.int64] = np.array([])

    if args.json == None:
        inJCoups_origin: npt.NDArray[np.float64] = copy.deepcopy(inJCoups)

        while (1):
            # Delete Too Small JCoups J = args.lw*(-0.3) ~ args.lw*(0.3) use matrix Filter
            # 1: keep and 0: neglect
            mat_filter_low_factor: npt.NDArray[np.int64] = np.zeros(
                (inSParams.size, inSParams.size), dtype=np.int64)
            mat_filter_low_factor = (np.abs(inJCoups_origin) > args.thr) * 1
            inJCoups = copy.deepcopy(inJCoups_origin)
            inJCoups *= mat_filter_low_factor
            mat_filter_ab_quartet: npt.NDArray[np.int64] = np.zeros(
                (inSParams.size, inSParams.size), dtype=np.int64)

            # np.fill_diagonal(mat_filter_ab_quartet, 0)
            import math
            for idx, x in enumerate(inSParams):
                for idy, y in enumerate(inSParams):
                    if idx == idy:
                        mat_filter_ab_quartet[idx][idy] = 0
                    else:
                        # if two chemical shift is very close , will perform AB quartet
                        # if x-y == 0 the Ratio_J_Hz will crash
                        # if the JCoups is negative number, will first priority to use normal QM cal.
                        # if not in AB Quartet and negative nubmers will use Multiplet (nmrsim)
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
                                if mat_filter_low_factor[idx][idy] == 0:
                                    pass
                                else:
                                    ic(idx, x, idy, y)
                                    print("Something wrong in your SParams !!!")
                                    ic()
                                    raise ValueError(
                                        idx + x + idy + y + " was not found or is a directory")

            # ic(mat_filter_ab_quartet)

            mat_filter_multi = mat_filter_low_factor - mat_filter_ab_quartet

            idx0_ab_connect: list[list[int | set[int]]] = []
            for idx, x in enumerate(mat_filter_ab_quartet):
                np_nonzero: npt.NDArray[np.int64] = (x*(idx+1)).nonzero()[0]
                group: set[int] = set(np_nonzero.tolist())
                # group = set(((x*(idx+1)).nonzero()[0]))
                group.add(idx)
                idx0_ab_connect.append([idx, group])

            idx0_ab_group_set = []
            for _, x in idx0_ab_connect:
                if len(x) == 1:                                 # type: ignore
                    idx0_ab_group_set.append(x)
                elif len(x) > 1:                                # type: ignore
                    group: set[int] = x
                    jump: bool = False
                    bond_penetration: int = 1
                    while jump == False:
                        jump = True
                        for idy, y in enumerate(group):
                            if (not group.issuperset(idx0_ab_connect[y][1])) and bond_penetration <= args.tb:
                                jump = False
                                group = group.union(idx0_ab_connect[y][1])
                            bond_penetration += 1
                    idx0_ab_group_set.append(group)
                else:
                    print("Something wrong in ab_group_set")
                    ic()
                    raise ValueError("idx0_ab_group_set is bugs !!!")

            # CH3 Equivalent is manually control by symmetry
            # So if chemical shift in AB quartet region need to move to multiplet
            list_Equivalent3: list[int] = []
            for idx, x in enumerate(inAnmr.nucinfo[1]):
                if x[1] == 3:
                    for idy, y in enumerate(in_idxAtoms):
                        if y == min(x[2]):
                            list_Equivalent3.append(idy)
            set_Equivalent3: set[int] = set(list_Equivalent3)

            # Equivalent3 is idx0 numbers
            for idx, x in enumerate(idx0_ab_group_set):
                set_move: set[int] = x.intersection(set_Equivalent3)
                if not len(set_move) == 0:
                    idx0_ab_group_set[idx] = set(x).difference(
                        set_move).union(set([idx]))
                    set_move = set_move.difference(set([idx]))
                for idy, y in enumerate(set_move):
                    mat_filter_multi[idx][y] = 1

            print(" ===== Processing =====")
            print(" threshold of JCoupling  : ", f'{args.thr:>3.5f}')
            print(" threshold of AB quartet : ", f'{args.thrab:>3.5f}')
            print("  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")

            max_len_AB: int = 0
            for idx, x in enumerate(idx0_ab_group_set):
                print(f'{(idx+1):>5d}{len(x):>5d}', {a+1 for a in x}, set(
                    a+1 for a in [y*idy for idy, y in enumerate(mat_filter_multi[idx].tolist()) if y != 0]).difference({a+1 for a in x}))
                if len(x) > max_len_AB:
                    max_len_AB = len(x)
            if (max_len_AB > args.mss):
                if args.auto:
                    args.thrab = args.thrab + 0.0025
                    args.thr = args.thr + 0.030
                else:
                    print(" Need to tidy the nspins of AB quartet and use cmd -thrab ")
                    exit(0)
            else:
                # print(" Use this parameter to calculate the Full Spectra")
                break
        #
        # AB quartet if more than two peaks of AB quartet, added closed peaks (not in AB quartet)
        print(" ===== Modification AB quartet =====")
        for idx, idx0_set_origin in enumerate(idx0_ab_group_set):
            idx0_list: list[int] = list(idx0_set_origin)
            # print(set(inSParams[list(idx0_set_origin)]))
            if len(set(list(inSParams[idx0_list]))) == 1:
                from censo_ext.Tools.spectra import find_nearest
                # ic(inSParams[list(idx0_set_origin)])
                _, Move_idx0 = find_nearest(list(inSParams[list(idx0_list)]),
                                            inSParams[idx0_list[0]])
                a = np.argwhere(
                    inSParams[:] == inSParams[int(Move_idx0)])
                idx0_ab_group_set[idx] = idx0_set_origin.union(
                    set(int(idx) for idx in a[0]))
                # print(set(list(inSParams[list(idx0_set_origin)])))
        for idx, idx0_set_origin in enumerate(idx0_ab_group_set):
            idx0_list = list(idx0_set_origin)
            print(f'{(idx+1):>5d}{len(idx0_list):>5d}', {a+1 for a in idx0_list}, set(
                a+1 for a in [idx0_set*idx for idx, idx0_set in enumerate(mat_filter_multi[idx].tolist()) if idx0_set != 0]).difference({a+1 for a in idx0_list}))
        print(" Use this parameter to calculate the Full Spectra")

    # Low level QM model
    # see https://nmrsim.readthedocs.io/en/latest/index.html
    #
    import censo_ext.Tools.qm as qm
    idx0_peaks_range: list[int] = []
    # ic(args.lw)

    if args.json == None:
        print("")
        print(" ===== Processing =====")
        print(" the group of calculate spectra :", len(idx0_ab_group_set))
        print("  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")
        Res_peaks: list[list[tuple[float, float]]] = []
        for idx, idx0_set_origin in enumerate(idx0_ab_group_set):
            idx0_list: list[int] = list(idx0_set_origin)
            print(f'{(idx+1):>5d}{len(idx0_list):>5d}', {a+1 for a in idx0_list}, set(
                a+1 for a in [idx0_set*idx for idx, idx0_set in enumerate(mat_filter_multi[idx].tolist()) if idx0_set != 0]).difference({a+1 for a in idx0_list}))

            v: npt.NDArray[np.float64] = inSParams[idx0_list]
            J: npt.NDArray[np.float64] = inJCoups[idx0_list].T[idx0_list]

            Res_qm_base: list[tuple[float, float]] = qm.qm_base(v=list(
                v), J=J, nIntergals=inHydrogen[idx0_list.index(idx)], idx0_nspins=idx0_list.index(idx), args=args)

            Res_qm_multiplet: list[tuple[float, float]] = []
            for _, z in enumerate(Res_qm_base):
                list_multiplicity: list[npt.NDArray[np.int64]] = list(set([x*idx for idx, x in enumerate(
                    mat_filter_multi[idx]) if x != 0]).difference({a for a in idx0_list}))
                inJ: list[tuple[float, int]] = []
                for ida, a in enumerate(list_multiplicity):
                    if np.fabs(inSParams[idx]-inSParams[a]) > 0.1:
                        inJ.append((inJCoups[idx][a], inHydrogen[a]))

                if len(inJ) >= 1:
                    tmp: npt.NDArray[np.float64] = np.array(
                        qm.qm_multiplet(z[0], nIntergals=1, J=inJ))
                    tmp.T[1] *= z[1]
                    Res_qm_multiplet += tmp.tolist()
                elif len(inJ) == 0:
                    Res_qm_multiplet = Res_qm_base
                else:
                    print("Something wrong")
                    ic()
                    exit(1)
            from nmrsim.math import normalize_peaklist
            Res_qm_multiplet: list[tuple[float, float]] = np.array(
                normalize_peaklist(Res_qm_multiplet, inHydrogen[idx])).tolist()

            if len(Res_peaks) == 0 and len(Res_qm_multiplet) == 0:
                pass
            else:
                Res_peaks.append(Res_qm_multiplet)

        import json
        import pickle
        with open(str(inAnmr.get_Directory()/Path("peaks.json")), "w") as final:
            json.dump(Res_peaks, final)

        idx0_peaks_range = [*range(len(Res_peaks))]

    else:
        import json
        import pickle
        with open(str(inAnmr.get_Directory()/Path("peaks.json")), "r") as json_file:
            Res_peaks = json.load(json_file)

        if args.json[0] == -1:
            idx0_peaks_range = [*range(len(Res_peaks))]
        else:
            idx0_peaks_range = args.json

    Final_peaks: list[tuple[float, float]] = []
    for idx, x in enumerate(Res_peaks):
        if idx in idx0_peaks_range:
            Final_peaks += x

    print("")
    print(" ===== Processing to plotting spectra =====")
    print(" Wait a minutes ...")
    args.out = str(inAnmr.get_Directory()/Path(args.out))

    if dpi != None and Active_range != None:
        np_dat: npt.NDArray[np.float64] = qm.print_plot(
            Final_peaks, dpi, nIntergals=1, args=args, Active_range=Active_range, hidden=not args.show)
    else:
        print("dpi and Active_range is wrong")
        ic()
        exit(0)

    if not args.bobyqa:
        print(" the spectra is saved to : ", args.out)

    print(" All done ...")
    return np_dat


if __name__ == "__main__":
    main()
