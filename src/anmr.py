#!/usr/bin/env python3
import os
import copy
import numpy as np
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
        default=False,
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


def main(args: argparse.Namespace = argparse.Namespace()) -> np.ndarray:

    if args == argparse.Namespace():
        args = cml("")
    if args.dir == None:
        Directory_str: Path = Path(".")
    else:
        Directory_str: Path = Path(args.dir)

    from Tools.anmrfile import ClassAnmr
    inAnmr: ClassAnmr = ClassAnmr(Directory_str)
    inAnmr.method_read_anmrrc()
    inAnmr.method_read_nucinfo()
    # ic(inAnmr.get_average_orcaSJ_Exist(), args.average)
    if inAnmr.get_average_orcaSJ_Exist() and args.average:
        inAnmr.method_load_average_orcaSJ(args=args)
    else:
        inAnmr.method_read_enso()
        inAnmr.method_read_folder_orcaSJ()
        inAnmr.method_update_equivalent_orcaSJ()
        inAnmr.method_average_orcaSJ()
        inAnmr.method_save_average_orcaSJ()

    inSParams: np.ndarray = np.array(
        list(inAnmr.Average_orcaSJ.orcaSParams.values()))*args.mf

    # idxinSParams:np.ndarray = np.array(list(inAnmr.Average_orcaSJ.orcaSParams.keys()))
    # np.savetxt("inSParams.out",np.array(np.stack([idxinSParams,inSParams]).T), fmt=' %4.0f   %10.6f')

    inJCoups: np.ndarray = np.array(inAnmr.Average_orcaSJ.orcaJCoups)
    inidxAtoms: dict[int, str] = inAnmr.Average_orcaSJ.idxAtoms
    dpi: int | None = None
    Active_range: int | None = None
    inHydrogen: list[int] = []
    in_xyzFileName: Path = Path("crest_conformers.xyz")
    # ic(inAnmr.Directory)
    for idx, x in enumerate(inAnmr.get_Anmr_Active()):
        if idx == 0:
            if x == 'C':
                from Tools.ml4nmr import read_mol_neighbors_bond_order
                mol, neighbors, bond_order = read_mol_neighbors_bond_order(
                    inAnmr.get_Directory() / in_xyzFileName)
                inHydrogen = [value+1 for value in bond_order.values()]
                Active_range = 200
                dpi = 500
                args.lw = 20
                if not args.thr:
                    args.thr = args.lw * 0.3
                inJCoups = np.zeros_like(inJCoups)
                # ic()
                # os._exit(0)
            elif x == 'H':
                inHydrogenDict: dict = {}
                for key, value in inAnmr.Average_orcaSJ.idxAtoms.items():
                    inHydrogenDict[key] = len(inAnmr.nucinfo[1][key-1][2])
                for x in inAnmr.get_list_idx1AcidAtomsNoShowRemoveH(
                        inAnmr.get_Directory()/in_xyzFileName):
                    if x in inHydrogenDict.keys():
                        del inHydrogenDict[x]
                        del inidxAtoms[x]
                # ic(inHydrogenDict, len(inHydrogenDict))
                inHydrogen = [value for value in inHydrogenDict.values()]
                Active_range = 10
                dpi = 10000
                args.lw = 1
                if not args.thr:
                    args.thr = args.lw * 0.3
            else:
                print("Other Active Nuclear element, waiting to build")
                ic()
                os._exit(0)
        elif idx >= 1:
            print("only for ONE Active Nuclear element, waiting to build")
            ic()
            os._exit(0)

    # ic(inHydrogen,len(inHydrogen))
    # ic(inidxAtoms,len(inidxAtoms))
    # ic(inSParams,len(inSParams))
    # ic(inJCoups,inJCoups.shape)

    idx0_ab_group_set: list = []
    mat_filter_multi: np.ndarray = np.array([])
    # ic(args.json)
    if args.json == None:

        inJCoups_origin: np.ndarray = copy.deepcopy(inJCoups)
        # ic(inJCoups_origin)
        while (1):
            # Delete Too Small JCoups J = args.lw*(-0.3) ~ args.lw*(0.3) use matrix Filter
            # 1: keep and 0: neglect
            mat_filter_low_factor: np.ndarray = np.zeros(
                (inSParams.size, inSParams.size), dtype=int)
            mat_filter_low_factor = (np.abs(inJCoups_origin) > args.thr) * 1
            inJCoups = copy.deepcopy(inJCoups_origin)
            inJCoups *= mat_filter_low_factor
            # np.savetxt(str(inAnmr.Directory/Path("mat_filter_low_factor.out")),
            #           mat_filter_low_factor, fmt="%d")
            # np.savetxt(str(inAnmr.Directory/Path("JJ-filter.out")),
            #           inJCoups, fmt="%6.2f")
            mat_filter_ab_quartet: np.ndarray = np.zeros(
                (inSParams.size, inSParams.size), dtype=int)

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
                                    print("something wrong !!!")
                                    ic()
                                    os._exit(0)

            mat_filter_multi = mat_filter_low_factor - mat_filter_ab_quartet

            # if mat_filter_multi.size != 0 and mat_filter_low_factor.size != 0 and mat_filter_ab_quartet.size != 0:
            #    np.savetxt(str(inAnmr.Directory/Path("mat_filter_low_factor.out")),
            #               mat_filter_low_factor, fmt="%d")
            #    np.savetxt(str(inAnmr.Directory/Path("mat_filter_ab_quartet.out")),
            #               mat_filter_ab_quartet, fmt="%d")
            #    np.savetxt(
            #        str(inAnmr.Directory/Path("mat_filter_multi.out")), mat_filter_multi, fmt="%d")
            # else:
            #    ic()
            #    os._exit(0)

            idx0_ab_connect: list = []
            for idx, x in enumerate(mat_filter_ab_quartet):
                group = set((x*(idx+1)).nonzero()[0])
                group.add(idx)
                idx0_ab_connect.append([idx, group])

            idx0_ab_group_set = []
            for x in idx0_ab_connect:
                if len(x[1]) == 1:
                    idx0_ab_group_set.append(x[1])
                elif len(x[1]) > 1:
                    group = x[1]
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
                    ic()
                    os._exit(0)
            # ic(idx0_ab_group_set)

            # CH3 Equivalent is manually control by symmetry
            # So if chemical shift in AB quartet region need to move to multiplet
            List_Equivalent3: list = []
            for idx, x in enumerate(inAnmr.nucinfo[1]):
                if x[1] == 3:
                    for idy, y in enumerate(inidxAtoms):
                        if y == min(x[2]):
                            List_Equivalent3.append(idy)
            Set_Equivalent3 = set(List_Equivalent3)

            # Set_Equivalent3 is idx0 numbers
            for idx, x in enumerate(idx0_ab_group_set):
                Set_move = x.intersection(Set_Equivalent3)
                if not len(Set_move) == 0:
                    idx0_ab_group_set[idx] = set(x).difference(
                        Set_move).union(set([idx]))
                    Set_move = Set_move.difference(set([idx]))
                for idy, y in enumerate(Set_move):
                    mat_filter_multi[idx][y] = 1

            # np.savetxt("mat_filter_ab_quartet.out", mat_filter_ab_quartet, fmt="%d")
            # np.savetxt("mat_filter_multi.out", mat_filter_multi, fmt="%d")

            print(" ===== Processing =====")
            print(" threshold of JCoupling  : ", f'{args.thr:>3.5f}')
            print(" threshold of AB quartet : ", f'{args.thrab:>3.5f}')
            print("  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")

            max_len_AB: int = 0
            for idx, x in enumerate(idx0_ab_group_set):
                print(f'{(idx+1):>5d}{len(x):>5d}', {a+1 for a in x}, set(
                    a+1 for a in [y*idy for idy, y in enumerate(mat_filter_multi[idx]) if y != 0]).difference({a+1 for a in x}))
                if len(x) > max_len_AB:
                    max_len_AB = len(x)
            if (max_len_AB > args.mss):
                if args.auto:
                    args.thrab = args.thrab + 0.0025
                    args.thr = args.thr + 0.030
                else:
                    print(" Need to tidy the nspins of AB quartet and use cmd -thrab ")
                    ic()
                    os._exit(0)
            else:
                print(" Use this parameter to calculate the Full Spectra")
                break

    # Low level QM model
    # see https://nmrsim.readthedocs.io/en/latest/index.html
    #
    import Tools.qm as qm
    idx0_peaks_range: list = []

    if args.json == None:
        print("")
        print(" ===== Processing =====")
        print(" the group of calculate spectra :", len(idx0_ab_group_set))
        print("  idx len(x) {x's AB quartet} {x's all - x's AB quartet} ")
        Result_peaks: list = []
        for idx, idx0_set_origin in enumerate(idx0_ab_group_set):
            idx0_set: list = list(idx0_set_origin)
            print(f'{(idx+1):>5d}{len(idx0_set):>5d}', {a+1 for a in idx0_set}, set(
                a+1 for a in [idx0_set*idx for idx, idx0_set in enumerate(mat_filter_multi[idx]) if idx0_set != 0]).difference({a+1 for a in idx0_set}))

            v: np.ndarray = inSParams[idx0_set]
            J: np.ndarray = inJCoups[idx0_set].T[idx0_set]
            # ic(v,J)
            Result_qm_base: list[np.ndarray] = qm.qm_base(v=list(
                v), J=J, nIntergals=inHydrogen[idx0_set.index(idx)], idx0_nspins=idx0_set.index(idx), args=args)
            # ic(idx, Result_qm_base)
            Result_qm_multiplet: list[np.ndarray] = []
            for idz, z in enumerate(Result_qm_base):
                list_multiplicity: list = list(set([x*idx for idx, x in enumerate(
                    mat_filter_multi[idx]) if x != 0]).difference({a for a in idx0_set}))
                inJ: list = []
                for ida, a in enumerate(list_multiplicity):
                    import math
                    if math.fabs(inSParams[idx]-inSParams[a]) > 0.1:
                        inJ.append((inJCoups[idx][a], inHydrogen[a]))

                if len(inJ) >= 1:
                    tmp: np.ndarray = np.array(
                        qm.qm_multiplet(z[0], nIntergals=1, J=inJ))
                    tmp.T[1] *= z[1]
                    Result_qm_multiplet += tmp.tolist()
                elif len(inJ) == 0:
                    Result_qm_multiplet = Result_qm_base
                else:
                    ic()
                    os._exit(0)
            from nmrsim.math import normalize_peaklist
            Result_qm_multiplet: list[np.ndarray] = np.array(
                normalize_peaklist(Result_qm_multiplet, inHydrogen[idx])).tolist()

            if len(Result_peaks) == 0 and len(Result_qm_multiplet) == 0:
                pass
            else:
                Result_peaks.append(Result_qm_multiplet)

        import json
        import pickle
        with open(str(inAnmr.get_Directory()/Path("peaks.json")), "w") as final:
            json.dump(Result_peaks, final)

        idx0_peaks_range: list = [*range(len(Result_peaks))]

        # with open("peaks_setting.json","wb") as final:
        #    pickle.dump(idx0_ab_group_set,final)

    else:
        import json
        import pickle
        with open(str(inAnmr.get_Directory()/Path("peaks.json")), "r") as json_file:
            Result_peaks = json.load(json_file)

        # with open("peaks_setting.json","rb") as pickle_file:
        #    idx0_ab_group_set = pickle.load(pickle_file)
        if args.json[0] == -1:
            idx0_peaks_range: list = [*range(len(Result_peaks))]
        else:
            idx0_peaks_range: list = args.json

    Final_Result_peaks: list[np.ndarray] = []
    for idx, x in enumerate(Result_peaks):
        if idx in idx0_peaks_range:
            Final_Result_peaks += x

    print("")
    print(" ===== Processing to plotting spectra =====")
    print(" Wait a minutes ...")
    args.out = str(inAnmr.get_Directory()/Path(args.out))

    if dpi != None and Active_range != None:
        np_dat: np.ndarray = qm.print_plot(
            Final_Result_peaks, dpi, nIntergals=1, args=args, Active_range=Active_range, hidden=not args.show)
        # np_dat: np.ndarray = qm.print_plot(
        #    Final_Result_peaks, dpi, nIntergals=1, args=args, Active_range=Active_range, hidden=not args.show)
    else:
        print("dpi and Active_range is wrong")
        ic()
        os._exit(0)

    if not args.bobyqa:
        print(" the spectra is saved to : ", args.out)

    print(" All done ...")
    return np_dat


if __name__ == "__main__":
    main()
