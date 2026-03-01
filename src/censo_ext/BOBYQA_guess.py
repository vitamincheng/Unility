#!/usr/bin/env python
import argparse
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import nmrglue as ng
from censo_ext.Tools.anmrfile import AD_Normal
from censo_ext.Tools.datfile import CensoDat, Peaks_npz, unit_conversion
from censo_ext.Tools.utility import AtomID, IsExist_bool, print_arguments

descr = """
_______________________________________________________________________________
| For generate the peak.npz for orcaS-BOBYQA.out or Integral of spectra  
| Usages    : BOBYQA_guess.py <geometry> [options]
| [options]
| File      : -i input dat/npz file [default 1r.npz]
| Auto      : --atuo Automated mode and read dat/npz file [default False]
| Basic     : --basic Only one time for threshold under automated mode [default False]
| Manual    : -m --manual Manual mode and read the peaks.npz [default False]
| Save      : --save To save peaks.npz [default False]
| Show      : -show --show Show spectra on screen [default false]
| threshold : -t -thr threshold of peaks [default 1.0]
| phase     : -p --phase phase of spectra (1 to -1) [default 1.0]
| Delete    : --delete Delete specific cID peaks
| Merge     : --merge Merge cID peaks to one peak
| Cut       : --cut Cut cID peak to two peaks by lowest point
|______________________________________________________________________________
"""

useit = """\
    End     Endremove    Startremove                 Start
    +               +    +                               +
    +---------------+----+-------------------------------+
    lower field                               higher field
                        delta /ppm
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
        default="1r.npz",
        help="Provide input one npz file [default 1r.npz]",
    )

    parser.add_argument(
        "--auto",
        dest="auto",
        action="store_true",
        help="Automatically integrate the peaks [default False]",
    )

    parser.add_argument(
        "--basic",
        dest="basic",
        action="store_true",
        help="Only one time for threshold under automated mode [default False]",
    )

    parser.add_argument(
        "--save",
        dest="save",
        action="store_true",
        help="Saved the automatically information of the peaks to input file [default False]",
    )

    parser.add_argument(
        "-show",
        "--show",
        dest="show",
        action="store_true",
        help="Show the spectra on screen [default False]",
    )

    parser.add_argument(
        "-t",
        "--thr",
        dest="thr",
        action="store",
        type=float,
        required=False,
        default=1.0,
        help="threshold of peaks [default 1.0]",
    )

    parser.add_argument(
        "-p",
        "--phase",
        dest="phase",
        action="store",
        type=float,
        required=False,
        default=1.0,
        help="phase of spectra (1 to -1) [default 1.0]",
    )

    parser.add_argument(
        "-m",
        "--manual",
        dest="manual",
        action="store_true",
        help="Manually integrate and use plot_1D_peaks.out file [default False]",
    )

    parser.add_argument(
        "--delete",
        dest="delete",
        action="store",
        type=int,
        nargs="+",
        default=False,
        help="Delete specific cID peaks",
    )

    parser.add_argument(
        "--merge",
        dest="merge",
        action="store",
        type=int,
        nargs="+",
        default=False,
        help="Merge cID peaks to one peak",
    )

    parser.add_argument(
        "--cut",
        dest="cut",
        action="store",
        type=int,
        default=False,
        help="Cut cID peak to two peaks by lowest point",
    )

    return parser.parse_args()


# global variable
peaks_fileName = "peaks.npz"


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    if args.auto and args.manual:
        print("  For both args.auto and args.manual only for one mode")
        print("  Exit and Close the program !!!")
        exit(0)

    if not IsExist_bool(args.file):
        return

    censo: CensoDat = CensoDat(args.file)
    in_Data = censo.get_Dat().T
    ppm: npt.NDArray[np.float64] = in_Data[0]
    intensit: npt.NDArray[np.float64] = in_Data[1]

    y_heighest: float = float(np.max(intensit))
    y_lowest: float = float(np.min(intensit))

    from censo_ext.Tools.spectra import numpy_thr_mean_3
    thres: float = numpy_thr_mean_3(intensit)*args.thr + y_heighest * 0.01
    thres_baseline: float = thres
    uc: unit_conversion = unit_conversion(ppm)

    args_start, args_end = uc.ppm_limits()

    # plot and indicate all peaks
    fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
    ax: Axes = fig.subplots()
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                        top=0.90, wspace=0.05, hspace=0.05)

    peaks: Peaks_npz = Peaks_npz(uc)

    ng_1r_peaks = ng.peakpick.pick(
        data=intensit, pthres=thres, algorithm="downward")

    # Automatically Integate the peaks
    if args.auto:

        thres, ng_1r_peaks = process_auto_mode(
            args, intensit, y_heighest, thres, thres_baseline, uc, peaks, ng_1r_peaks)

    # Integrate the peaks if manually fixed the peaks.npz file
    if args.manual:
        peaks.method_read_file()
        print("  ========== Before ==========")
        peaks.method_print()

        # Delete peaks if manually read the file
        if args.delete:
            peaks.method_delete_cID(args.delete)

        # Merge peaks if manually read the file
        if args.merge:
            peaks.method_merge_cID(args.merge)

        # Cut peak if manually read the file
        if args.cut:
            peaks.method_cut_cID(args.cut, intensit)

        if args.delete or args.cut or args.merge:
            print("  ========== After ==========")
            peaks.method_print()

    # Draw the integral lines and cID of peaks
    # Plot the integration lines, limits and cID of peaks
    if args.auto or args.manual:
        draw_Integral(intensit, ax, peaks)

    # add markers for peak positions. It is only for preview.
    if not args.auto and not args.manual:
        draw_preview(intensit, uc, args_start, args_end, ax, ng_1r_peaks)

    # draw the threshold line and text and for adjust threshold for next time
    if args.auto:
        draw_threshold(thres, args_start, args_end, ax)

    if args.show:
        # draw the spectra
        draw_spectra(intensit, uc)

        # draw the x axis
        draw_axis(args, y_heighest, y_lowest, args_start, args_end, fig, ax)
        plt.show()

    if args.save:
        peaks.method_save()


def draw_axis(args: argparse.Namespace, y_heighest: float, y_lowest: float, start: float, end: float, fig: Figure, ax: Axes) -> None:
    plt.xlim(end, start)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", which="both", bottom=True,
                   top=False, labelbottom=True, labelsize=12)
    ax.tick_params(axis="y", which="both", left=False,
                   right=False, labelleft=False)
    ax.get_yaxis().set_visible(False)

    # If phase is -1, it will adjust the y axis
    if y_lowest*(-1) < y_heighest*0.2:
        plt.ylim(-0.05*y_heighest, 1.10*y_heighest)
    else:
        plt.ylim(1.10*y_lowest, 1.10*y_heighest)
    fig.suptitle(args.file, fontsize=12, y=0.98)
    fig.text(0.5, 0.04, "$\\delta$ / ppm", ha="center", fontsize=12)


def draw_spectra(intensit: npt.NDArray[np.float64], uc: unit_conversion) -> None:
    plt.plot(uc.ppm_scale(), intensit, 'b', linewidth=1)


def draw_threshold(thres: float, start: float, end: float, ax: Axes) -> None:
    plt.hlines(thres, end, start, linestyles="--")  # type: ignore # nopep8
    ax.text(start, thres*1.02, f"thr = {thres:>10.3f}",
            ha="center", va="center")


def draw_preview(intensit: npt.NDArray[np.float64], uc: unit_conversion, args_start: float, args_end: float, ax: Axes, ng_1r_peaks) -> None:
    idx_cID: float = 0

    for idx_peaks, cID, LW, VOL in ng_1r_peaks:
        if idx_cID < cID:
            idx_cID = cID
        else:
            break
        height: float = intensit[int(idx_peaks)]
        ppm_peak: float = uc.ppm(idx_peaks)

        min: int = uc.index(args_start)
        max: int = uc.index(args_end)
        if ppm_peak < max and ppm_peak > min:
            ax.scatter(ppm_peak, height, marker="o", color="r", s=100, alpha=0.5)  # type: ignore # nopep8
            ax.text(ppm_peak, height*1.05, str(cID),
                    ha="center", va="center")


def draw_Integral(intensit: npt.NDArray[np.float64], ax: Axes, peaks: Peaks_npz) -> None:
    Data: list[tuple[int, npt.NDArray[np.float64],
                     npt.NDArray[np.float64]]] = peaks.method_integrate(intensit)
    for cID, peak_int, peak_scale in Data:
        ax.plot(peak_scale, peak_int.cumsum() /
                100./4 + peak_int.max()*0.8, 'g-')
        ax.text(peak_scale[0], 0.5 * peak_int.sum() / 100./4 + peak_int.max()*0.8, str(cID),
                fontsize=8)


def process_auto_mode(args: argparse.Namespace, intensit, y_heighest, thres, thres_baseline, uc, peaks, ng_1r_peaks):

    AD_normal: AD_Normal = AD_Normal()
    if not AD_normal.Exist():
        print("Only use --basic for one parameter under Automated mode")
        args.basic = True
        nGroups = 0
    else:
        AD_normal.method_load_files()
        if isinstance(AD_normal.ChemicalShifts, dict):
            idx1_orcaS: list[int] = list(
                map(int, AD_normal.ChemicalShifts.keys()))
        else:
            print("  The format of orcaS.out is not dict !!!")
            print("  Exit and Close the program !!!")
            exit(0)
        from censo_ext.Tools.anmrfile import Anmr
        inAnmr: Anmr = Anmr()
        inAnmr.method_read_nucinfo()
        ChemEqvs: dict[AtomID, list[AtomID]] = {key: value for key,
                                                value in inAnmr.NeighborChemEqvs.items() if key in idx1_orcaS}
        Groups: list[list[int]] = list(
            sorted(value) for value in ChemEqvs.values())
        unique_group: list = []
        for item in Groups:
            if item not in unique_group:
                unique_group.append(item)
        nGroups: int = len(unique_group)

    peaks_list: list[tuple[int, float, float, float]] = []
    last_peaks: int = 0
    while True:

        peaks_list = extract_peak(intensit, uc, ng_1r_peaks)

        print("threshold : ", thres)
        print("Excepted  : ", nGroups)
        print("Real Num  : ", len(peaks_list))

        merge_overlap_peaks(peaks_list)

        for x in peaks_list:
            print(x)
        print("")
        peaks.method_load_Data(peaks_list)

        if args.basic is True or nGroups == len(peaks_list) or len(peaks_list) < last_peaks:
            break
        elif len(peaks_list) > last_peaks:
            last_peaks = len(peaks_list)
        else:  # find the smallest of len(peak_list)
            if thres > y_heighest*0.7:
                break
            else:
                thres += thres_baseline
                ng_1r_peaks = ng.peakpick.pick(
                    data=intensit, pthres=thres, algorithm="downward")

    print("  ========== Automated Data ==========")
    print("Numbers of peaks : ", len(peaks))
    peaks.method_print()
    return thres, ng_1r_peaks


def merge_overlap_peaks(peak_list) -> None:
    ppm_end: list[np.float64] = np.array(peak_list).T[2].tolist()
    ppm_start: list[np.float64] = np.array(peak_list).T[1].tolist()
    ppm_end.pop(0)
    ppm_end.append(999)  # type: ignore
    ppm_args: npt.NDArray[np.intp] = np.argwhere(
        np.array(ppm_end)-np.array(ppm_start) < 0)
    for x in (ppm_args+1):
        index: int = x[0]
        ppm_center: float = (peak_list[index-1]
                             [1] + peak_list[index][2])/2
        new_cID, start, end, Area = peak_list[index-1]
        peak_list[index-1] = (new_cID, ppm_center, end, Area)
        new_cID, start, end, Area = peak_list[index]
        peak_list[index] = (new_cID, start, ppm_center, Area)


def extract_peak(intensit, uc, ng_1r_peaks) -> list[tuple[int, float, float, float]]:
    peak_list: list[tuple[int, float, float, float]] = []
    sorted_cID_peaks: npt.NDArray = np.sort(
        ng_1r_peaks, order='cID')
    new_cID: list[int] = []
    for cID in sorted_cID_peaks['cID']:
        args_cID: npt.NDArray[np.intp] = (
            np.argwhere(sorted_cID_peaks['cID'] == cID))
        r_Axis: int = int(
            sorted_cID_peaks[args_cID.min()]['X_AXIS'])
        l_Axis: int = int(
            sorted_cID_peaks[args_cID.max()]['X_AXIS'])
        r_LW: float = float(
            sorted_cID_peaks[args_cID.min()]['X_LW'])
        l_LW: float = float(
            sorted_cID_peaks[args_cID.max()]['X_LW'])

        if l_LW <= 1.0:
            l_LW = 1
        if r_LW <= 1.0:
            r_LW = 1
        l_LW_thr: float = 120/l_LW
        r_LW_thr: float = 120/r_LW

        l_peak: float = uc.ppm(l_Axis)+(l_LW/3000)*l_LW_thr
        r_peak: float = uc.ppm(r_Axis)-(r_LW/3000)*r_LW_thr

        min: int = uc.index(l_peak)
        max: int = uc.index(r_peak)
        if min > max:
            min, max = max, min

            # extract the peak
        peak: npt.NDArray[np.float64] = intensit[min:max + 1]

        if cID not in new_cID:
            new_cID.append(cID)
            peak_list.append(
                (int(cID), l_peak, r_peak, float(peak.sum())))

    return peak_list


if __name__ == "__main__":
    main()
