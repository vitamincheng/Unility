#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
from icecream import ic
import matplotlib.pyplot as plt
import nmrglue as ng
from censo_ext.Tools.utility import IsExist_bool
descr = """
________________________________________________________________________________
|                                                      vitamin.cheng@gmail.com
| For generate orcaS.BOBYQA or Intergal of spectra  
| Usages   : BOBYQA_guess.py <geometry> [options]
| [options]
| File     : -i input npz file [default 1r.npz]
| Auto     : --atuo Automated [default False]
| Save     : --save To save peaks.npz [default False]
| Start    : -start Start point of chemical shift [default from data] 
| End      : -end End point of chemical shift [default from data]
| threshold: -t -thr threshold of peaks [default 1.0]
| phase    : -p --phase phase of spectra (1 to -1) [default 1.0]
| Manual   : -m --manual Manual mode [default False]
| Delete   : --delete Delete specific cID peaks
| Merge    : --merge Merge cID peaks to one peak
| Cut      : --cut Cut cID peak to two peaks by lowest point
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
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="1r.npz",
        help="Provide one input npz file [default 1r.npz]",
    )

    parser.add_argument(
        "--auto",
        dest="auto",
        action="store_true",
        help="Automatically integrate the peaks [default False]",
    )

    parser.add_argument(
        "--save",
        dest="save",
        action="store_true",
        help="Saved the automatically information of the peaks to input file [default False]",
    )

    parser.add_argument(
        "-start",
        dest="start",
        action="store",
        type=float,
        required=False,
        default=None,
        help="start point of chemical shift [default from data]",
    )

    parser.add_argument(
        "-end",
        dest="end",
        action="store",
        type=float,
        required=False,
        default=None,
        help="end point of chemical shift [default from data]",
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

    args: argparse.Namespace = parser.parse_args()
    return args


# global variable
peaks_fileName = "peaks.npz"


class unit_conversion():

    def __init__(self, in_ppm: npt.NDArray[np.float64]) -> None:
        self.in_ppm: npt.NDArray[np.float64] = in_ppm
        self.args_ppm: dict = {idx: ppm for idx,
                               ppm in enumerate(in_ppm)}
        self.ppm_args: dict = {ppm: idx for idx,
                               ppm in enumerate(in_ppm)}

    def index(self, ppm) -> int:
        from censo_ext.Tools.spectra import find_nearest
        value, index = find_nearest(list(self.in_ppm), ppm)
        return index

    def ppm(self, index) -> float:
        return float(self.args_ppm[index])

    def ppm_scale(self) -> npt.NDArray[np.float64]:
        return self.in_ppm

    def ppm_limits(self) -> tuple[float, float]:
        return float(self.in_ppm.min()), float(self.in_ppm.max())


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()

    if IsExist_bool(args.file):

        in_Data = np.load(args.file)["arr_0"]
        a, *b = in_Data.shape
        if a == 1:
            in_Data = in_Data[0]
        ppm: npt.NDArray[np.float64] = in_Data.T[0]
        intensit: npt.NDArray[np.float64] = in_Data.T[1]

        y_heighest: np.float64 = np.max(intensit)
        y_lowest: np.float64 = np.min(intensit)

        from censo_ext.Tools.spectra import numpy_thr_mean_3
        threshold: float = numpy_thr_mean_3(intensit)*args.thr
        threshold += y_heighest * 0.01

        uc: unit_conversion = unit_conversion(ppm)

        if args.start is None or args.end is None:
            args.start, args.end = uc.ppm_limits()
        elif args.start and args.end:
            if args.start > args.end:
                args.start, args.end = args.end, args.start
        else:
            print("  Args.start or args.end have wroing !!!")
            print("  Exit and Close to the program !!!")
            exit(0)

        # plot and indicate all peaks
        fig = plt.figure(figsize=(11.7, 8.3), dpi=100)
        ax = fig.subplots()
        fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                            top=0.90, wspace=0.05, hspace=0.05)

        np_peaks: npt.NDArray = np.array(
            [], dtype=[('cID', 'i8'), ('Start', 'f8'), ('End', 'f8'), ('Area', 'f8')])

        ng_peaks = ng.peakpick.pick(
            data=intensit, pthres=threshold, algorithm="downward")

        # Automatically Intergate the peaks
        if args.auto:
            peak_list: list = []
            sorted_cID_peaks: npt.NDArray = np.sort(ng_peaks, order='cID')
            new_cID: list[int] = []
            for cID in sorted_cID_peaks['cID']:
                args_cID: npt.NDArray[np.intp] = (
                    np.argwhere(sorted_cID_peaks['cID'] == cID))
                r_Axis: float = float(
                    sorted_cID_peaks[args_cID.min()]['X_AXIS'])
                l_Axis: float = float(
                    sorted_cID_peaks[args_cID.max()]['X_AXIS'])
                r_LW: float = float(sorted_cID_peaks[args_cID.min()]['X_LW'])
                l_LW: float = float(sorted_cID_peaks[args_cID.max()]['X_LW'])

                if l_LW <= 1.0:
                    l_LW = 1
                if r_LW <= 1.0:
                    r_LW = 1
                l_LW_thr: float = 120/l_LW
                r_LW_thr: float = 120/r_LW

                l_peak: float = uc.ppm(l_Axis)+(l_LW/10000)*l_LW_thr
                r_peak: float = uc.ppm(r_Axis)-(r_LW/10000)*r_LW_thr

                min: int = uc.index(l_peak)
                max: int = uc.index(r_peak)
                if min > max:
                    min, max = max, min

                # extract the peak
                peak: npt.NDArray[np.float64] = intensit[min:max + 1]

                if cID not in new_cID:
                    new_cID.append(cID)
                    peak_list.append((cID, l_peak, r_peak, peak.sum()))

            np_peaks = np.array(
                peak_list, dtype=[('cID', 'i8'), ('Start', 'f8'), ('End', 'f8'), ('Area', 'f8')])

        # Intergrate the peaks if manually fixed the peaks.npz file
        if args.manual:
            in_Peaks = np.load(peaks_fileName)["arr_0"]
            ic(in_Peaks)
            a, *b = in_Peaks.shape
            if a == 2:
                np_peaks = in_Peaks[1][0]
            else:
                np_peaks = in_Peaks[0]

            idx: int = 0
            for cID, start, end, _ in np_peaks:  # type: ignore
                min: int = uc.index(start)
                max: int = uc.index(end)
                if min > max:
                    min, max = max, min

                # extract the peak
                peak = intensit[min:max + 1]
                peak_scale = uc.ppm_scale()[min:max + 1]

                np_peaks['Area'][idx] = peak.sum()  # type: ignore
                idx += 1

            # Delete or Merge peaks if manually read the file
            if args.delete:
                for x in args.delete:
                    if x in np_peaks['cID']:
                        np_peaks = np_peaks[np_peaks['cID'] != x]
                    else:
                        print("your delete element is wrong cID")
                        exit(0)

            if args.merge:
                min_cID: int = args.merge.min()
                start, end = -99999, 99999
                for x in sorted(args.merge):
                    if x in np_peaks['cID']:
                        args_x: npt.NDArray[np.intp] = np.argwhere(
                            np_peaks['cID'] == x)
                        if start < np_peaks[args_x[0]]['Start']:
                            start: float = np_peaks[args_x[0]
                                                    ]['Start'].max().astype(float)
                        if end > np_peaks[args_x[0]]['End']:
                            end: float = np_peaks[args_x[0]
                                                  ]['End'].min().astype(float)
                    else:
                        print("  Your merge element is wrong cID")
                        exit(0)

                # assign the new cID of data
                args_x = np.argwhere(np_peaks['cID'] == min_cID)
                min: int = uc.index(start)
                max: int = uc.index(end)
                Total_intensit = 0
                for x in args.merge:
                    Total_intensit += np.sum(
                        np_peaks[np_peaks['cID'] == x]['Area'])

                # Remove unnecessary entry
                for x in args.merge:
                    if (x in np_peaks['cID']):
                        np_peaks = np_peaks[np_peaks['cID'] != x]

                np_peaks = np.insert(
                    np_peaks, args_x[0], (min_cID, start, end, Total_intensit))

            if args.cut:

                if args.cut in np_peaks['cID']:
                    # use ng.peakpick.pick from y_heighest 0.99 to down to two different peaks
                    args_x = np.argwhere(np_peaks['cID'] == args.cut)
                    l_peaks: float = np_peaks[args_x][0]['Start'][0].astype(
                        float)
                    r_peaks: float = np_peaks[args_x][0]['End'][0].astype(
                        float)
                    min: int = uc.index(l_peaks)
                    max: int = uc.index(r_peaks)
                    if min > max:
                        min, max = max, min
                    cut_peaks = ng.peakpick.pick(
                        data=intensit[min:max+1], pthres=threshold, algorithm="downward")
                    sorted_cut_peaks = np.sort(cut_peaks, order='VOL')

                    end = int(sorted_cut_peaks['X_AXIS'][-1] + min)
                    start = int(sorted_cut_peaks['X_AXIS'][-2] + min)
                    cut_argmin: np.intp = np.argmin(intensit[start:end+1])
                    cut_center: float = uc.ppm(start+cut_argmin)

                    # remove the old entry and add two additional entry
                    np_peaks = np_peaks[np_peaks['cID'] != args.cut]

                    np_peaks = np.insert(
                        np_peaks, args_x[0], (args.cut, cut_center, r_peaks, intensit[min:start+cut_argmin].sum()))
                    np_peaks = np.insert(
                        np_peaks, args_x[0], (np_peaks['cID'].max()+1, l_peaks, cut_center, intensit[start+cut_argmin:max].sum()))
                else:
                    print("your merge element is wrong cID")
                    ic()
                    exit(1)

        # Draw the intergral lines
        if args.auto or args.manual:
            for cID, start, end, _ in np_peaks:  # type: ignore
                min: int = uc.index(start)
                max: int = uc.index(end)
                if min > max:
                    min, max = max, min

                # extract the peak
                peak = intensit[min:max + 1]
                peak_scale: npt.NDArray[np.float64] = uc.ppm_scale()[
                    min:max + 1]

                # plot the integration lines, limits and name of peaks
                ax.plot(peak_scale, peak.cumsum() /
                        100./4 + peak.max()*0.8, 'g-')
                # ic(cID, peak.sum())
                # ax.plot(peak_scale, [0] * len(peak_scale), 'r-')
                ax.text(peak_scale[0], 0.5 * peak.sum() / 100./4 + peak.max()*0.8, cID,
                        fontsize=8)

        # add markers for peak positions. It is only for preview
        if not args.auto and not args.manual:
            idx_cID: float = 0
            peak_list: list = []
            for idx_peaks, cID, LW, VOL in ng_peaks:
                if idx_cID < cID:
                    idx_cID = cID
                else:
                    break
                height: float = intensit[int(idx_peaks)]
                ppm_peak: float = uc.ppm(idx_peaks)

                min: int = uc.index(args.start)
                max: int = uc.index(args.end)
                if ppm_peak < max and ppm_peak > min:
                    ax.scatter(ppm_peak, height, marker="o", color="r", s=100, alpha=0.5)  # type: ignore # nopep8
                    ax.text(ppm_peak, height*1.05, str(cID),
                            ha="center", va="center")

        # draw the spectra
        plt.plot(uc.ppm_scale(), intensit, 'b', linewidth=1)

        # draw the threshold line and text
        if not args.manual:
            plt.hlines(threshold, args.end, args.start, linestyles="--")  # type: ignore # nopep8
            ax.text(args.start, threshold*1.02, f"thr = {threshold:>10.3f}",
                    ha="center", va="center")

        # draw the x axis
        plt.xlim(args.end, args.start)
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

        plt.show()

        if args.save:
            from censo_ext.Tools.utility import save_simulation_spectra_file_npz
            save_simulation_spectra_file_npz(peaks_fileName, np_peaks)

        # if len(np_peaks) > 0:
        #    np.savetxt("test.out", np_peaks, fmt="%10i %18f %18f %18f")


if __name__ == "__main__":
    main()
