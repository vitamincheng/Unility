#! /usr/bin/env python
from pathlib import Path
from icecream import ic
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import nmrglue as ng
import matplotlib.pyplot as plt
from nmrglue.fileio.fileiobase import unit_conversion
import numpy as np
import numpy.typing as npt
import argparse
import sys
from censo_ext.Tools.spectra import numpy_thr_mean_3
from censo_ext.Tools.utility import delete_all_files, print_arguments


descr = """
________________________________________________________________________________
| For Plot 1D sepctra in experiments using nmrglue module
| Usage: plot_1D_DEPT.py <geometry> [options]
| [Options]
| Input    : -i the pdata path(under 1r folder) [required]
|          : -start start point of chemical shift [default from data]
|          : -end   end point of chemical shift [default from data]
| Save     : --save saved the report of carbon [default false]
| Hidden   : --hidden show the plot [default False]
|______________________________________________________________________________
"""
useit = """
    End     Endremove    Startremove                 Start
    +               +    +                               +
    +---------------+----+-------------------------------+
    lower field                               higher field
                        delta /ppm
    """


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=f"{descr}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
        add_help=True
    )

    parser.add_argument(
        "-d",
        "--dir",
        dest="dir",
        action="store",
        type=str,
        required=False,
        help="Provide the the parent's path of your 13C / DEPT_90 / DEPT_135 including the pdata",
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
        "--save",
        dest="save",
        action="store_true",
        help="Saved the report of carbon [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


# global variable
pipe_fid_filename = ".1d_pipe.fid"
peaks_fileName = "plot_1D_DEPT.peaks"


class diagram:
    def __init__(self, path: dict[Path, Path], args: argparse.Namespace) -> None:
        self._fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
        self._ax_list: list[Axes] = self._fig.subplots(
            3, 1, sharex=True)  # type: ignore
        self._fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                                  top=0.90, wspace=0.05, hspace=0.05)
        self._data: dict[Path, npt.NDArray] = {}
        self._uc_1h: dict[Path, unit_conversion] = {}
        self._path: dict[Path, Path] = path
        self._limits: tuple[float, float] = args.end, args.start
        self._thr: dict[Path, float] = {}
        self._title = "Normal Mode"
        self._key = ""
        self._ax_selected = ""
        self.load_data()
        self.draw_curve()
        self.draw_x_axis(full=True)
        self.draw_title()

    def load_data(self):
        print(self._path)
        for x, value in self._path.items():
            self.load_singlet_data(x, 1)
        self.load_singlet_data(Path("DEPT_135"), -1)

    def load_singlet_data(self, x, phase):
        dic, data = ng.bruker.read_pdata(str(self._path[Path(x)]))
        udic = ng.bruker.guess_udic(dic, data)
        C = ng.convert.converter()
        C.from_bruker(dic, data, udic)
        ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
        dic, data = ng.pipe.read(pipe_fid_filename)
        data = data.real*phase  # type: ignore
        uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)
        if phase == 1:
            self._data[x] = data
            self._uc_1h[x] = uc_1h
        elif phase == -1:
            self._data[Path("DEPT_135_down")] = data
            self._uc_1h[Path("DEPT_135_down")] = uc_1h

    def draw_curve(self):
        self.draw_singlet_curve(
            self._data[Path('13C')], self._uc_1h[Path('13C')], self._ax_list[2])
        self.draw_singlet_curve(
            self._data[Path('DEPT_90')], self._uc_1h[Path('DEPT_90')], self._ax_list[0])
        self.draw_singlet_curve(
            self._data[Path('DEPT_135')], self._uc_1h[Path('DEPT_135')], self._ax_list[1])

    def draw_singlet_curve(self, data, uc_1h, ax):
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="y", which="both", left=False,
                       right=False, labelleft=False)
        ax.tick_params(axis="x", which="both", bottom=False,
                       top=False, labelbottom=False)
        ax.plot(uc_1h.ppm_scale(), data, 'b', linewidth=1)

    def draw_x_axis(self, full: bool = False):
        end, start = self._uc_1h[Path('13C')].ppm_limits()
        xlimits = self._ax_list[2].get_xlim()
        if isinstance(self._limits[0], (int, float)) and isinstance(self._limits[1], (int, float)):
            start = self._limits[1]
            end = self._limits[0]
        if full:
            self._ax_list[2].set_xlim(end, start)
        else:
            self._ax_list[2].set_xlim(xlimits[0], xlimits[1])
        self._ax_list[2].spines["bottom"].set_visible(True)
        self._ax_list[2].tick_params(axis="x", which="both", bottom=True,
                                     top=False, labelbottom=True, labelsize=12)
        self._fig.text(0.5, 0.04, "$\\delta$ / ppm",
                       ha="center", fontsize=12)

    def clear_local_axes(self) -> None:
        for ax in self._ax_list:
            _x, _y = ax.get_xlim(), ax.get_ylim()
            ax.clear()
            ax.set_xlim(_x)
            ax.set_ylim(_y)

    def on_key_press(self, event):
        """Callback function for key press events."""
        if event.key == 'q':
            print("Quitting the application.")
            plt.close(event.canvas.figure)

        elif event.key == 'h':
            self._show_help()

        elif event.key == 's':
            self._save_file()

        elif event.key == 't':
            self._set_threshold()

        elif event.key == 'd':
            self._display_carbon()

    def _save_file(self):
        pass

    def _set_threshold(self):
        self._title = "Set Threshold"
        self._key = "t"
        self.draw_title()
        self._fig.canvas.draw_idle()

    def _display_carbon(self):
        pass

    def on_button_release(self, event):
        from matplotlib.backend_bases import MouseButton
        toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore
        if self._key == "t" and toolbar_mode == "zoom rect":
            self._ax_list[0].set_navigate_mode("ZOOM")
            self._ax_list[1].set_navigate_mode("ZOOM")
            self._ax_list[2].set_navigate_mode("ZOOM")

        elif event.inaxes in self._ax_list and self._key == "t" and event.button == MouseButton.LEFT:
            self._button_xy = event.xdata, event.ydata
            # self._ax_selected = event.inaxes
            self._ax_selected = None
            self._fig.canvas.draw_idle()

    def on_button_press(self, event):
        from matplotlib.backend_bases import MouseButton
        toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore
        if self._key == "t" and toolbar_mode == "zoom rect":
            self._title = ""
            self._status_int = []
            # self.draw_title()
            self._ax_list[0].set_navigate_mode("ZOOM")
            self._ax_list[1].set_navigate_mode("ZOOM")
            self._ax_list[2].set_navigate_mode("ZOOM")
            self._fig.canvas.draw_idle()
        elif event.inaxes in self._ax_list and self._key == "t" and event.button == MouseButton.LEFT:
            self._button_xy = event.xdata, event.ydata
            self._ax_selected = event.inaxes
            start, end = event.inaxes.get_xlim()

            self.clear_local_axes()
            self.draw_curve()
            self.draw_title()
            self.draw_x_axis()
            event.inaxes.hlines(self._button_xy[1], end, start, colors='k', linestyles="dashed", linewidth=1)  # type: ignore # nopep8
            self._fig.canvas.draw_idle()

    def on_mouse_motion(self, event):
        from matplotlib.backend_bases import MouseButton

        if event.inaxes is self._ax_selected and event.button is MouseButton.LEFT:
            toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore

            if self._key == "t" and toolbar_mode == "zoom rect":
                self._title_bottom = ""
                self.draw_title()
                self._key = ""
                self._ax_list[0].set_navigate_mode("ZOOM")
                self._ax_list[1].set_navigate_mode("ZOOM")
                self._ax_list[2].set_navigate_mode("ZOOM")
                self._fig.canvas.draw_idle()
            if self._key == "t" and toolbar_mode != "zoom rect":

                self.clear_local_axes()
                self.draw_curve()
                self.draw_title()
                self.draw_x_axis()
                # self.draw_status()
                # self.draw_scatter_numbers()

                if event.inaxes in self._ax_list and self._key == "t" and event.inaxes == self._ax_selected:
                    start, end = event.inaxes.get_xlim()
                    # event.inaxes.hlines(self._button_xy[1], end, start, colors='k', linestyles="dashed", linewidth=1)  # type: ignore # nopep8
                    event.inaxes.hlines(event.ydata, end, start, colors='k', linestyles="dashed", linewidth=1)  # type: ignore # nopep8

                self._fig.canvas.draw_idle()

    def draw_title(self) -> None:
        """Draw title on the plot."""
        self._ax_list[0].set_title(self._title, loc="left")

    def _show_help(self):
        pass
        """Display help information."""
        print("Help mode ")
        self._status_int = []
        self._title = ("Press 'q' to Quit, 's' to Save file.\n"
                       "Press 't' to set the threshold.\n"
                       "Press 'd' to display the Carbon.\n")
        self.draw_title()
        self._fig.canvas.draw_idle()

    def connect(self):
        """Connect all event handlers."""
        self._cID_key = self._fig.canvas.mpl_connect(
            'key_press_event', self.on_key_press)
        self._cID_button_press = self._fig.canvas.mpl_connect(
            'button_press_event', self.on_button_press)
        self._cID_button_motion = self._fig.canvas.mpl_connect(
            'motion_notify_event', self.on_mouse_motion)
        self._cID_button_release = self._fig.canvas.mpl_connect(
            'button_release_event', self.on_button_release)

    def disconnect(self):
        self._fig.canvas.mpl_disconnect(self._cID_key)
        self._fig.canvas.mpl_disconnect(self._cID_button_press)
        self._fig.canvas.mpl_disconnect(self._cID_button_release)
        self._fig.canvas.mpl_disconnect(self._cID_button_motion)


def Channel(args, path, ax: Axes, phase: float = 1.0) -> dict:
    dic, data = ng.bruker.read_pdata(str(path))
    udic = ng.bruker.guess_udic(dic, data)
    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(pipe_fid_filename)
    data = data.real*phase  # type: ignore
    uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)

    # from censo_ext.Tools.spectra import numpy_thr_mean_3
    # threshold: float = numpy_thr_mean_3(data.astype(np.float64))*thr
    if phase == 1:
        # ax.hlines(threshold, args.end, args.start,
        #          linestyles="dashdot", linewidth=0.5)
        # if "DEPT_135" in str(path):
        #    threshold: float = numpy_thr_mean_3(
        #        data.astype(np.float64))*thr_ch3_180*(-1)
        # ax.hlines(threshold, args.end, args.start,
        #          linestyles="dashdot", linewidth=0.5)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", which="both", bottom=False,
                       top=False, labelbottom=False)
        ax.tick_params(axis="y", which="both", left=False,
                       right=False, labelleft=False)
        ax.set_xlim(args.end, args.start)
        ax.plot(uc_1h.ppm_scale(), data, 'b', linewidth=1)

    # end ---------+--------- start
    # args.end                args.start
    # ppm_1h_0                ppm_1h_1

    # ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    # ppm = np.linspace(ppm_1h_0, ppm_1h_1, data.shape[0])

    # if isinstance(args.start, (int, float)) and isinstance(args.end, (int, float)):
    #    pass
    # else:
    #    args.end, args.start = uc_1h.ppm_limits()

    # from censo_ext.Tools.spectra import numpy_thr_mean_3
    # threshold: float = 0
    # if isinstance(thr, float):
    #    threshold: float = numpy_thr_mean_3(data.astype(np.float64))*thr
    # detect all peaks with a threshold
    # from scipy.signal import find_peaks
    # y_heighest = max(data)
    # threshold += y_heighest * 0.01
    # peaks, _ = find_peaks(data, height=threshold, width=1)

    # add markers for peak positions
    Result: dict[int, float] = dict()

    # for n, peak in enumerate(peaks):
    #    ppm: float = uc_1h.ppm(peak)
    #    Result[n+1] = ppm

    return Result  # type: ignore


def Compare_two_dict(CH1: dict, CH2: dict, StAtoms: dict, Label: int) -> None:
    from censo_ext.Tools.spectra import find_nearest
    for x in CH2.values():
        nearest_peak, idx0 = find_nearest(list(CH1.values()), x)
        # ic(nearest_peak-x)
        if (nearest_peak-x) < 0.05:
            # if StAtoms[idx0+1] == -1:
            StAtoms[idx0+1] = Label
        else:
            print("some peaks is more than 0.02 ppm")
            print("  Exit and Close the program !!!")
            ic()
            exit(0)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    # read in the Bruker data
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    plt.rcParams['keymap.save'].remove('s')
    plt.rcParams['keymap.fullscreen'].remove('f')
    plt.rcParams['keymap.back'].remove('c')
    plt.rcParams['keymap.grid'].remove('g')
    plt.rcParams['toolbar'] = 'toolbar2'
    plt.ion()

    path: dict[Path, Path] = {}

    ch1: Path = Path('13C')
    ch2: Path = Path('DEPT_90')
    ch3: Path = Path('DEPT_135')

    # thr: dict[Path, float] = {ch1: 2.0, ch2: 20.0, ch3: 2.0}
    # thr_ch3_180: float = 2.0

    import platform
    _system = platform.system()
    if _system == "Linux":
        directory: Path = Path(
            "/home/vitamin/Simulation/38.Ergocalciferol(Vitamin_D2)/00.Spectra/bmse000510/nmr/set01")
    elif _system == "Darwin":
        directory: Path = Path(
            "/Users/chengwen-cheng/Desktop/Simulation/bmse000510/nmr/set01")
    else:
        print("  Only for ubuntu or Darwin system ...")
        print("  Exit and Close the program !!!")
        exit(0)

    path[ch1] = directory / ch1 / Path("pdata/1")
    path[ch2] = directory / ch2 / Path("pdata/1")
    path[ch3] = directory / ch3 / Path("pdata/1")

    # plot and indicate all peaks
    digrams: diagram = diagram(path, args)
    digrams.connect()
    plt.ioff()
    plt.show()
    delete_all_files(pipe_fid_filename)
    exit(0)
    # channel: Path = ch1
    # Result_ch1 = Channel(args, path=path[channel], ax=ax[2])

    StAtoms: dict[int, int] = {key: -1 for key, value in Result_ch1.items()}

    channel: Path = ch2
    Result_ch2 = Channel(args, path=path[channel], ax=ax[0])

    channel: Path = ch3
    Result_ch3 = Channel(args, path=path[channel], ax=ax[1])

    channel: Path = ch3
    Result_ch3_180 = Channel(args, path=path[channel], ax=ax[2], phase=-1.0)  # nopep8

    Compare_two_dict(Result_ch1, Result_ch3, StAtoms, Label=3)
    Compare_two_dict(Result_ch1, Result_ch2, StAtoms, Label=1)

    print("DEPT90             ppm")
    for idx1, ppm in Result_ch2.items():
        print(f"{idx1:6d} {ppm:>15.5f}")

    print("DEPT135(up)        ppm")
    for idx1, ppm in Result_ch3.items():
        print(f"{idx1:6d} {ppm:>15.5f}")

    Compare_two_dict(Result_ch1, Result_ch3_180, StAtoms, Label=2)

    print("DEPT135(down)      ppm")
    #
    for idx1, ppm in Result_ch3_180.items():
        print(f"{idx1:6d} {ppm:>15.5f}")

    for key, value in StAtoms.items():
        if value == -1:
            StAtoms[key] = 0

    # channel: Path = ch1
    # dic, data = ng.bruker.read_pdata(str(path[channel]))
    # udic = ng.bruker.guess_udic(dic, data)
    # C = ng.convert.converter()
    # C.from_bruker(dic, data, udic)
    # pipe_fid_filename = ".1d_pipe.fid"
    # ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
    # dic, data = ng.pipe.read(pipe_fid_filename)
    # data = data.real  # type: ignore
    # uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)

    # end ---------+--------- start
    # args.end                args.start
    # ppm_1h_0                ppm_1h_1

    # ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    # ppm: npt.NDArray[np.float64] = np.linspace(
    #    ppm_1h_0, ppm_1h_1, data.shape[0])

    # if isinstance(args.start, (int, float)) and isinstance(args.end, (int, float)):
    #    pass
    # else:
    #    args.end, args.start = uc_1h.ppm_limits()

    # threshold: float = numpy_thr_mean_3(
    #    data.astype(np.float64))*thr[channel]

    # detect all peaks with a threshold
    # from scipy.signal import find_peaks
    # y_heighest = max(data)
    # y_lowest = min(data)
    # threshold += float(y_heighest) * 0.01
    # peaks, _ = find_peaks(data, height=threshold, width=1)

    # print the final data
    # print("#   ID             ppm    nHydrogens")
    # for n, peak in enumerate(peaks):
    #    height = data[int(peak)]
    #    ppm = uc_1h.ppm(peak)
    #    print(f"{n+1:6d} {ppm:>15.5f}        {StAtoms[n+1]:>3d}")

    # save to file
    if args.save:
        with open(peaks_fileName, "w") as f:
            sys.stdout = f
            print("#   ID             ppm    nHydrogens")
            # for n, peak in enumerate(peaks):
            #    height = data[int(peak)]
            #    ppm = uc_1h.ppm(peak)
            #    print(f"{n+1:6d} {ppm:>15.5f}        {StAtoms[n+1]:>3d}")
            sys.stdout = sys.__stdout__

    # add markers for peak positions
    # for n, peak in enumerate(peaks):
    #    height = data[int(peak)]
    #    ppm = uc_1h.ppm(peak)
    #    if ppm < args.end and ppm > args.start:
    #        if StAtoms[n+1] == 0:
    #            ax[2].text(ppm.tolist(), height*1.20,
    #                       str("C"), ha="center", va="center", rotation=90)
    #        elif StAtoms[n+1] == 1:
    #            ax[2].text(ppm.tolist(), height*1.20,
    #                       "CH", ha="center", va="center", rotation=90)
    #        else:
    #            ax[2].text(ppm.tolist(), height*1.20,
    #                       "CH"+rf'$_{str(StAtoms[n+1])}$', ha="center", va="center", rotation=90)

    # fig.suptitle(args.path, fontsize=12, y=0.98)
    fig.text(0.5, 0.04, "$\\delta$ / ppm",
             ha="center", fontsize=12)
    ax[2].tick_params(axis="x", which="both", bottom=True,
                      top=False, labelbottom=True, labelsize=12)

    digrams.connect()
    plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()
