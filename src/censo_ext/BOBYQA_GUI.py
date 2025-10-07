#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import nmrglue as ng
from censo_ext.Tools.anmrfile import AD_Normal
from censo_ext.Tools.datfile import CensoDat, Peaks_npz, unit_conversion
from censo_ext.Tools.utility import IsExist_bool, print_arguments

descr = """
_______________________________________________________________________________
| For generate the peak.npz for orcaS-BOBYQA.out or Integral of spectra
| Usages    : BOBYQA_GUI.py <geometry> [options]
| [options]
| File      : -i input dat/npz file [default 1r.npz]
| Auto      : --atuo Automated mode and read dat/npz file [default False]
| Basic     : --basic Only one time for threshold under automated mode [default False]
| Manual    : -m --manual Manual mode and read the peaks.npz [default False]
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


class diagram():

    # Fig = plt.figure(figsize=(11.7, 8.3), dpi=100)
    # Ax = fig.subplots()
    # Fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
    #                    top=0.90, wspace=0.05, hspace=0.05)
    def __init__(self, plt, fileName, peaks_npz, intensit, uc, thres, ng_1r_peaks) -> None:
        # self._key: str = 'e'
        self._plt = plt
        self._start, self._end = uc.ppm_limits()
        self._fileName = fileName
        self._peaks_npz = peaks_npz
        self._intensit = intensit
        self._button_x = None
        self._button_y = None
        self._uc: unit_conversion = uc
        self._thres = thres
        self._ng_1r_peaks = ng_1r_peaks
        self._fig = plt.figure(figsize=(11.7, 8.3), dpi=100)
        self._ax = self._fig.subplots()
        self._fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                                  top=0.90, wspace=0.05, hspace=0.05)
        self._ax.set_title("Edit mode.\nPress 'h' to help. ", loc="left")
        self._click_button: list[int] = []

    def on_key_press(self, event):
        """Callback function for key press events."""
        if event.key == 'q':
            print("Quitting the application.")
            plt.close(event.canvas.figure)

        elif event.key == 'enter':
            xmin, xmax, ymin, ymax = self._plt.axis()
            self._click_button = list(map(int, set(self._click_button)))
            self._click_button.sort()
            print(self._click_button)
            if self._key == "d":
                self._peaks_npz.method_delete_cID(self._click_button)
            elif self._key == "m":
                self._peaks_npz.method_merge_cID(self._click_button)
            elif self._key == "c":
                if len(self._click_button) == 1:
                    self._peaks_npz.method_cut_cID(
                        self._click_button[0], self._intensit)
                else:
                    print("  length of click button under cut mode more than 1 ")
                    print("  Exit and Close the program !!!")
                    exit(0)

            self._ax.clear()
            self._plt.xlim(xmin, xmax)
            self._plt.ylim(ymin, ymax)
            self.draw_threshold()
            self.draw_curve()
            self.draw_integral()
            self._click_button = []
            self._fig.canvas.draw_idle()

        elif event.key == 'escape':
            print("Edit mode")
            self._ax.set_title("Edit mode.\nPress 'h' to help. ", loc="left")
            self._key = 'e'
            self._fig.canvas.draw_idle()
        elif event.key == 'h':
            print("Help mode")
            self._ax.set_title(
                "Press 'q' to quit, 's' to save file.\nPress 'm' to merge / 'd' to delete / 'c' to cut mode\nPress 'Esc' to Edit mode, 'Enter' to Executive mode", loc="left")
            self._fig.canvas.draw_idle()
        elif event.key == 's':
            self._ax.set_title("Save to peaks.npz file", loc="left")
            self._peaks_npz.method_save()
            self._fig.canvas.draw_idle()
        elif event.key == 'm':
            print("Merge mode : ", end="")
            self._click_button = []
            self._ax.set_title("Merge mode", loc="left")
            self._key = 'm'
            self._fig.canvas.draw_idle()
        elif event.key == 'd':
            print("delete mode: ", end="")
            self._click_button = []
            self._ax.set_title("Delete mode", loc="left")
            self._key = 'd'
            self._fig.canvas.draw_idle()
        elif event.key == 'c':
            print("cut mode : ", end="")
            self._click_button = []
            self._ax.set_title("Cut mode", loc="left")
            self._key = 'c'
            self._fig.canvas.draw_idle()
        elif event.key == 'f':
            self._ax.clear()
            self.draw_x_axis()
            self.draw_threshold()
            self.draw_curve()
            self.draw_integral()
            self._fig.canvas.draw_idle()

    def on_button_release(self, event):
        if event.inaxes == self._ax and self._button_x is not None and self._button_y is not None:
            release_x = event.xdata
            release_y = event.ydata
            distance = np.sqrt((release_x-self._button_x) **
                               2 + (release_y-self._button_y)**2)
            tolerance = 0.01
            cID: int
            if distance < tolerance:
                Result = self._peaks_npz.method_ppm2cID(release_x)
                if Result is None:
                    return None
                else:
                    cID = int(Result)
            else:
                return None
            self._button_x = None
            self._button_y = None

            if self._key == 'c' and len(self._click_button) == 0:
                self._click_button.append(cID)
            elif self._key == 'c' and len(self._click_button) == 1:
                pass
            else:
                self._click_button.append(cID)

            self._ax.set_title(f"{self._click_button}",
                               loc="right", fontsize=10)

    def on_button_press(self, event):
        if event.inaxes == self._ax:
            if event.button == 1:
                self._button_x = event.xdata
                self._button_y = event.ydata

    def connect(self):
        self._cid_key = self._fig.canvas.mpl_connect(
            'key_press_event', self.on_key_press)
        self._cid_button = self._fig.canvas.mpl_connect(
            'button_press_event', self.on_button_press)
        self._cid_button_release = self._fig.canvas.mpl_connect(
            'button_release_event', self.on_button_release)

    def disconnect(self):
        self._fig.canvas.mpl_disconnect(self._cid_key)
        self._fig.canvas.mpl_disconnect(self._cid_button)
        self._fig.canvas.mpl_disconnect(self._cid_button_release)

    def draw_integral(self):
        Data = self._peaks_npz.method_integrate(self._intensit)
        for cID, peak_int, peak_scale in Data:
            self._ax.plot(peak_scale, peak_int.cumsum() /
                          100./5 + peak_int.max()*0.8, 'g-')
            self._ax.text(peak_scale[0], 0.5 * peak_int.sum() / 100./4 + peak_int.max()*0.8, cID,
                          fontsize=8)
        a = self._peaks_npz.get_cIDs_center_peaks()
        for ppm in a[1]:
            index = self._uc.index(ppm)
            height = self._intensit[index]
            self._ax.scatter(ppm, height, marker="o", color="r", s=30, alpha=0.5)  # type: ignore # nopep8

    def draw_preivew(self):

        idx_cID: float = 0
        for idx_peaks, cID, _, _ in self._ng_1r_peaks:
            if idx_cID < cID:
                idx_cID = cID
            else:
                break
            height: float = self._intensit[int(idx_peaks)]
            ppm_peak: float = self._uc.ppm(idx_peaks)

            args_start, args_end = self._uc.ppm_limits()
            min: int = self._uc.index(args_start)
            max: int = self._uc.index(args_end)
            if ppm_peak < max and ppm_peak > min:

                self._ax.scatter(ppm_peak, height, marker="o", color="r", s=100, alpha=0.5)  # type: ignore # nopep8
                self._ax.text(ppm_peak, height*1.05, str(cID),
                              ha="center", va="center")

    def draw_threshold(self):
        self._plt.hlines(self._thres, self._end, self._start, linestyles="--")  # type: ignore # nopep8
        self._ax.text(self._start, self._thres*1.02, f"thr = {self._thres:>10.3f}",
                      ha="center", va="center")

    def draw_curve(self):
        plt.plot(self._uc.ppm_scale(), self._intensit, 'b', linewidth=1)

    def draw_x_axis(self):

        y_heighest: float = float(np.max(self._intensit))
        y_lowest: float = float(np.min(self._intensit))
        plt.xlim(self._end, self._start)
        self._ax.spines["right"].set_visible(False)
        self._ax.spines["top"].set_visible(False)
        self._ax.spines["left"].set_visible(False)
        self._ax.tick_params(axis="x", which="both", bottom=True,
                             top=False, labelbottom=True, labelsize=12)
        self._ax.tick_params(axis="y", which="both", left=False,
                             right=False, labelleft=False)
        self._ax.get_yaxis().set_visible(False)

        # If phase is -1, it will adjust the y axis
        if y_lowest*(-1) < y_heighest*0.2:
            plt.ylim(-0.05*y_heighest, 1.10*y_heighest)
        else:
            plt.ylim(1.10*y_lowest, 1.10*y_heighest)
        self._fig.suptitle(self._fileName, fontsize=12, y=0.98)
        self._fig.text(0.5, 0.04, "$\\delta$ / ppm", ha="center", fontsize=12)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    plt.rcParams['keymap.save'].remove('s')
    plt.rcParams['keymap.fullscreen'].remove('f')
    plt.rcParams['keymap.back'].remove('c')
    plt.ion()

    if args.auto and args.manual:
        print("  For both args.auto and args.manual only for one mode")
        print("  Exit and Close the program !!!")
        exit(0)

    if IsExist_bool(args.file):

        censo: CensoDat = CensoDat(args.file)
        in_Data = censo.get_Dat().T
        ppm: npt.NDArray[np.float64] = in_Data[0]
        intensit: npt.NDArray[np.float64] = in_Data[1]

        y_heighest: float = float(np.max(intensit))
        from censo_ext.Tools.spectra import numpy_thr_mean_3
        thres: float = numpy_thr_mean_3(intensit)*args.thr
        thres += y_heighest * 0.01
        thres_baseline: float = thres
        uc: unit_conversion = unit_conversion(ppm)

        peaks_npz: Peaks_npz = Peaks_npz(uc)
        ng_1r_peaks = ng.peakpick.pick(
            data=intensit, pthres=thres, algorithm="downward")
        diagrams: diagram = diagram(plt, args.file,
                                    peaks_npz, intensit, uc, thres, ng_1r_peaks)

        # Automatically Integate the peaks
        if args.auto:

            AD_normal: AD_Normal = AD_Normal()
            if (AD_normal.Exist()):
                AD_normal.method_load_files()
                if isinstance(AD_normal.SParams, dict):
                    idx1_orcaS: list[int] = list(
                        map(int, AD_normal.SParams.keys()))
                else:
                    print("  The format of orcaS.out is not dict !!!")
                    print("  Exit and Close the program !!!")
                    exit(0)
                from censo_ext.Tools.anmrfile import Anmr
                inAnmr: Anmr = Anmr()
                inAnmr.method_read_nucinfo()
                ChemEqvs: dict[int, list[int]] = {key: value for key,
                                                  value in inAnmr.NeighborChemEqvs.items() if key in idx1_orcaS}
                Groups: list[list[int]] = list(
                    sorted(value) for value in ChemEqvs.values())
                unique_group: list = []
                for item in Groups:
                    if item not in unique_group:
                        unique_group.append(item)
                nGroups: int = len(unique_group)
            else:
                print("Only use --basic for one parameter under Automated mode")
                args.basic = True
                nGroups = 0

            peak_list: list = []
            last_peaks: int = 0
            while (1):
                peak_list = []
                sorted_cID_peaks: npt.NDArray = np.sort(
                    ng_1r_peaks, order='cID')
                new_cID: list[int] = []
                for cID in sorted_cID_peaks['cID']:
                    args_cID: npt.NDArray[np.intp] = (
                        np.argwhere(sorted_cID_peaks['cID'] == cID))
                    r_Axis: float = float(
                        sorted_cID_peaks[args_cID.min()]['X_AXIS'])
                    l_Axis: float = float(
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
                        peak_list.append(
                            (int(cID), l_peak, r_peak, float(peak.sum())))

                print("threshold : ", thres)
                print("Excepted  : ", nGroups)
                print("Real Num  : ", len(peak_list))
                print("     cID        Start          End             Area")
                for x in peak_list:
                    print(f"{x[0]:8d} {x[1]:12.4f} {x[2]:12.4f} {x[3]:16.4e}")
                print("")
                peaks_npz.method_load_Data(peak_list)

                if args.basic is True:
                    break
                if nGroups == len(peak_list):
                    break
                elif len(peak_list) < last_peaks:
                    break
                elif len(peak_list) > last_peaks:
                    last_peaks = len(peak_list)
                else:  # find the smallest of len(peak_list)
                    if thres > y_heighest*0.7:
                        break
                    else:
                        thres += thres_baseline
                        ng_1r_peaks = ng.peakpick.pick(
                            data=intensit, pthres=thres, algorithm="downward")

            print("  ========== Automated Data ==========")
            print("Numbers of peaks : ", len(peaks_npz))
            peaks_npz.method_print()

        # Integrate the peaks if manually fixed the peaks.npz file
        if args.manual:
            peaks_npz.method_read_file()
            print("  ========== Before ==========")
            peaks_npz.method_print()

        # Draw the intergral lines and cID of peaks
        # Plot the integration lines, limits and cID of peaks
        if args.auto or args.manual:
            diagrams.draw_integral()

        # add markers for peak positions. It is only for preview.
        if not args.auto and not args.manual:
            diagrams.draw_preivew()

        # draw the threshold line and text and for adjust threshold for next time
        if args.auto:
            diagrams.draw_threshold()

        diagrams.draw_curve()
        diagrams.draw_x_axis()
        diagrams.connect()

        plt.ioff()
        plt.show()


if __name__ == "__main__":
    main()
