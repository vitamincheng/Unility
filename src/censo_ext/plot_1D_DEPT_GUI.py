#! /usr/bin/env python
from enum import Enum
from pathlib import Path
from icecream import ic
from matplotlib.axes import Axes
from matplotlib.backend_bases import KeyEvent, MouseEvent
from matplotlib.figure import Figure
import nmrglue as ng
import matplotlib.pyplot as plt
from nmrglue.fileio.fileiobase import unit_conversion
import numpy as np
import numpy.typing as npt
import argparse
from censo_ext.Tools.utility import delete_all_files, print_arguments


descr = """
________________________________________________________________________________
| For Plot 1D sepctra in experiments using nmrglue module
| Usage: plot_1D_DEPT.py <geometry> [options]
| [Options]
| Input    : -i the pdata path(under 1r folder) [required]
|          : -start start point of chemical shift [default from data]
|          : -end   end point of chemical shift [default from data]
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
        # description=f"{descr}",
        # formatter_class=argparse.RawDescriptionHelpFormatter,
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

    return parser.parse_args()


# global variable
pipe_fid_filename = ".1d_pipe.fid"
peaks_fileName = "plot_1D_DEPT.peaks"


class DEPT(Enum):
    _13C = Path("13C")
    _90 = Path("DEPT_90")
    _135 = Path("DEPT_135")
    _135_down = Path("DEPT_135_down")


class diagram:
    def __init__(self, directory: Path, args: argparse.Namespace) -> None:
        self._fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
        self._ax_list: list[Axes] = self._fig.subplots(
            3, 1, sharex=True)  # type: ignore
        self._ax_dict: dict[DEPT, Axes] = {
            DEPT._13C: self._ax_list[2], DEPT._90: self._ax_list[0], DEPT._135: self._ax_list[1]}
        del self._ax_list
        self._fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                                  top=0.90, wspace=0.05, hspace=0.05)
        self._data: dict[DEPT, npt.NDArray] = {}
        self._uc_1h: dict[DEPT, unit_conversion] = {}
        self._directory: Path = directory
        self._limits: tuple[float, float] = args.end, args.start
        self._thr: dict[DEPT, float] = {
            DEPT._13C: 0.0, DEPT._90: 0.0, DEPT._135: 0.0, DEPT._135_down: 0.0}
        # index, and ppm
        self._Result: dict[DEPT, dict[int, float]] = {}
        self._title = "Normal Mode"
        self._key: str = ""
        self._ax_selected: Axes
        # index: numbers of Hyddrogen bonding to Atoms
        self._nHydrogens_Atom: dict[int, list[int]] = {}
        self._peaks: dict[int, int] = {}  # index of peaks in Data
        self.load_data()
        self.draw_curve()
        self.draw_x_axis(full=True)
        self.draw_title()

    def load_data(self) -> None:
        self.load_singlet_data(DEPT._13C, 1)
        self.load_singlet_data(DEPT._90, 1)
        self.load_singlet_data(DEPT._135, 1)
        self.load_singlet_data(DEPT._135, -1)

    def load_singlet_data(self, x: DEPT, phase: float) -> None:
        pdata = Path("pdata/1")
        dic, data = ng.bruker.read_pdata(
            str(self._directory / x.value / pdata))
        udic = ng.bruker.guess_udic(dic, data)
        C = ng.convert.converter()
        C.from_bruker(dic, data, udic)
        ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
        dic, data = ng.pipe.read(pipe_fid_filename)
        data: npt.NDArray = data.real*phase  # type: ignore
        uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)
        if phase == 1:
            self._data[x] = data
            self._uc_1h[x] = uc_1h
        elif phase == -1:
            self._data[DEPT._135_down] = data
            self._uc_1h[DEPT._135_down] = uc_1h

    def draw_curve(self) -> None:
        for x in [DEPT._13C, DEPT._90, DEPT._135]:
            self.draw_singlet_curve(
                self._data[x], self._uc_1h[x], self._ax_dict[x])

    def draw_singlet_curve(self, data, uc_1h, ax) -> None:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="y", which="both", left=False,
                       right=False, labelleft=False)
        ax.tick_params(axis="x", which="both", bottom=False,
                       top=False, labelbottom=False)
        ax.plot(uc_1h.ppm_scale(), data, 'b', linewidth=1)

    def draw_x_axis(self, full: bool = False) -> None:
        end, start = self._uc_1h[DEPT._13C].ppm_limits()
        xlimits: tuple[float, float] = self._ax_dict[DEPT._13C].get_xlim()
        if isinstance(self._limits[0], (int, float)) and isinstance(self._limits[1], (int, float)):
            start: float | int = self._limits[1]
            end: float | int = self._limits[0]
        if full:
            plt.xlim(end, start)
            self._limits = end, start
        else:
            self._ax_dict[DEPT._13C].set_xlim(xlimits[0], xlimits[1])
        self._ax_dict[DEPT._13C].spines["bottom"].set_visible(True)
        self._ax_dict[DEPT._13C].tick_params(axis="x", which="both", bottom=True,
                                             top=False, labelbottom=True, labelsize=12)
        self._fig.text(0.5, 0.04, "$\\delta$ / ppm",
                       ha="center", fontsize=12)

    def draw_status(self) -> None:
        if self._key == "t":
            # display the Total numbers
            for x in [DEPT._13C, DEPT._90]:
                self._ax_dict[x].set_title(
                    f"Threshold : {self._thr[x]:7.2f}\n", loc="right", x=1.05, y=0, fontsize=8)
            self._ax_dict[DEPT._135].set_title(
                f"Threshold : {self._thr[DEPT._135]:7.2f}\n\
                Threshold : {self._thr[DEPT._135_down]:7.2f}\n", loc="right", x=1.05, y=0, fontsize=8)

        if self._key == "d":
            # display the Total numbers
            for x in [DEPT._13C, DEPT._90]:
                self._ax_dict[x].set_title(
                    f"Total Numbers : {len(self._Result[x])}\n", loc="right", x=1.05, y=0, fontsize=8)
            self._ax_dict[DEPT._135].set_title(
                f"Total Numbers : {len(self._Result[DEPT._135])}\n\
                Total Numbers : {len(self._Result[DEPT._135_down])}\n", loc="right", x=1.05, y=0, fontsize=8)

            # display the total Hydrogens
            xmin, xmax = self._ax_dict[DEPT._13C].get_xlim()
            ymin, ymax = self._ax_dict[DEPT._13C].get_ylim()
            y: float = (ymax-ymin)*0.05
            self._ax_dict[DEPT._13C].text(
                xmin, y, f"Total Hydrogens : {sum([x[0] for x in self._nHydrogens_Atom.values()])}", ha="center", va="center", rotation=0)

    def clear_local_axes(self) -> None:
        for ax in self._ax_dict.values():
            _x, _y = ax.get_xlim(), ax.get_ylim()
            ax.clear()
            ax.set_xlim(_x)
            ax.set_ylim(_y)

    def on_key_press(self, event: KeyEvent) -> None:
        """Callback function for key press events."""
        if event.key == 'q':
            print("Quitting the application.")
            plt.close(event.canvas.figure)

        elif event.key == 'escape':
            self._reset_to_edit_mode()

        elif event.key == 'h':
            self._show_help()

        elif event.key == 's':
            self._save_file()

        elif event.key == 'f':
            self._redraw_full()

        elif event.key == 't':
            self._set_threshold()

        elif event.key == 'd':
            self._display_carbon()

    def _redraw_full(self) -> None:
        xmin, xmax, ymin, ymax = plt.axis()
        for x in self._ax_dict.values():
            x.clear()
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        self.draw_curve()
        self.draw_title()
        self.draw_x_axis(full=True)
        self._fig.canvas.draw_idle()

    def _reset_to_edit_mode(self) -> None:
        print("Normal mode")
        self._key = 'escape'

        self.clear_local_axes()
        self._title = "Normal mode.\nPress 'h' to help. "
        self.draw_curve()
        self.draw_title()
        self.draw_x_axis()
        self._fig.canvas.draw_idle()

    def _save_file(self) -> None:
        out: list | npt.NDArray = []

        print("#   ID             ppm        Hydrogens     manual")
        for key, value in self._Result[DEPT._13C].items():
            print(
                f"{key:6d} {float(value):>15.5f} {self._nHydrogens_Atom[key][0]:>12d} {self._nHydrogens_Atom[key][1]:>12d}")
            out.append([float(key), float(value),
                       int(self._nHydrogens_Atom[key][0]), int(self._nHydrogens_Atom[key][1])])
        out = np.array(out)
        np.savetxt(peaks_fileName, out, fmt='%6d %15.5f %12d %12d',
                   header="ID             ppm       nHydrogens     manual", comments="  # ")
        print(f" Saved the file to {peaks_fileName}")

    def _set_threshold(self) -> None:
        self._title = "Set Threshold"
        self._key = "t"
        self.draw_title()
        self._fig.canvas.draw_idle()

    def _display_carbon(self) -> None:
        self._title = "Display Carbon"
        self._key = "d"

        self.clear_local_axes()
        self.draw_curve()
        self.draw_x_axis()
        self.draw_title()
        self.cal_singlet_ch(DEPT._13C)
        self.cal_singlet_ch(DEPT._90)
        self.cal_singlet_ch(DEPT._135)
        self.cal_singlet_ch(DEPT._135_down)
        self.cal_st()
        self.draw_carbon_number()
        self.draw_status()
        self._fig.canvas.draw_idle()

    def draw_carbon_number(self) -> None:

        for idx1, peak in self._peaks.items():
            height: float = self._data[DEPT._13C][int(peak)]
            ppm: float = self._uc_1h[DEPT._13C].ppm(peak)
            end, start = self._limits
            if ppm < end and ppm > start:
                AX: Axes = self._ax_dict[DEPT._13C]
                if self._nHydrogens_Atom[idx1][1] == 0:
                    _COLOR = "k"
                else:
                    _COLOR = "r"

                if self._nHydrogens_Atom[idx1][0] == 0:
                    AX.text(ppm, height*1.20,
                            str("C"), ha="center", va="center", rotation=90, color=_COLOR)
                elif self._nHydrogens_Atom[idx1][0] == 1:
                    AX.text(ppm, height*1.20,
                            "CH", ha="center", va="center", rotation=90, color=_COLOR)
                else:
                    AX.text(ppm, height*1.20,
                            "CH"+rf'$_{str(self._nHydrogens_Atom[idx1][0])}$', ha="center", va="center", rotation=90, color=_COLOR)

    def cal_st(self) -> None:
        self._nHydrogens_Atom = {
            key: [-1, 0] for key, _ in self._Result[DEPT._13C].items()}
        self.Compare_two_dict(DEPT._13C, DEPT._135, Label=3)
        self.Compare_two_dict(DEPT._13C, DEPT._90, Label=1)
        self.Compare_two_dict(DEPT._13C, DEPT._135_down, Label=2)

        for key, value in self._nHydrogens_Atom.items():
            if value[0] == -1:
                self._nHydrogens_Atom[key][0] = 0

    def Compare_two_dict(self, CH1: DEPT, CH2: DEPT, Label: int) -> None:
        from censo_ext.Tools.spectra import find_nearest
        for x in self._Result[CH2].values():
            nearest_peak, idx0 = find_nearest(
                list(self._Result[CH1].values()), x)
            if (nearest_peak-x) < 0.05:
                self._nHydrogens_Atom[idx0+1][0] = Label
            else:
                print("some peaks is more than 0.02 ppm")
                print("  Exit and Close the program !!!")
                ic()
                exit(0)

    def cal_singlet_ch(self, ch: DEPT) -> None:
        thr: float = self._thr[ch]
        data = self._data[ch]
        uc_1h: unit_conversion = self._uc_1h[ch]
        from scipy.signal import find_peaks
        peaks, _ = find_peaks(data, height=thr, width=1)

        # add markers for peak positions
        Result: dict[int, float] = dict()

        for n, peak in enumerate(peaks, 1):
            ppm: float = uc_1h.ppm(peak)
            Result[n] = ppm
        self._Result[ch] = Result
        if ch == DEPT._13C:
            for n, peak in enumerate(peaks, 1):
                self._peaks[n] = peak

    def on_button_release(self, event: MouseEvent) -> None:
        from matplotlib.backend_bases import MouseButton
        toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore
        if self._key == "t" and toolbar_mode == "zoom rect":
            for x in self._ax_dict.values():
                x.set_navigate_mode("ZOOM")

        elif event.ydata is not None and event.inaxes in self._ax_dict.values() and self._key == "t" and event.button == MouseButton.LEFT:
            if event.inaxes == self._ax_dict[DEPT._90]:
                self._thr[DEPT._90] = event.ydata
            elif event.inaxes == self._ax_dict[DEPT._135] and event.ydata >= 0:
                self._thr[DEPT._135] = event.ydata
            elif event.inaxes == self._ax_dict[DEPT._135] and event.ydata < 0:
                self._thr[DEPT._135_down] = -event.ydata
            elif event.inaxes == self._ax_dict[DEPT._13C]:
                self._thr[DEPT._13C] = event.ydata
            self.draw_status()
            self._fig.canvas.draw_idle()

    def on_button_press(self, event: MouseEvent) -> None:
        from matplotlib.backend_bases import MouseButton
        toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore
        if self._key == "t" and toolbar_mode == "zoom rect":
            self._title = ""
            for x in self._ax_dict.values():
                x.set_navigate_mode("ZOOM")
            self._fig.canvas.draw_idle()
        elif event.inaxes is not None and event.inaxes in self._ax_dict.values() and self._key == "t" and event.button == MouseButton.LEFT:

            if event.xdata is not None and event.ydata is not None:
                self._button_xy: tuple[float, float] = event.xdata, event.ydata

            self._ax_selected: Axes = event.inaxes
            start, end = event.inaxes.get_xlim()

            self.clear_local_axes()
            self.draw_curve()
            self.draw_title()
            self.draw_x_axis()
            self.draw_status()
            event.inaxes.hlines(self._button_xy[1], end, start, colors='k', linestyles="dashed", linewidth=1)  # type: ignore # nopep8
            self._fig.canvas.draw_idle()

        elif event.xdata is not None and event.inaxes is self._ax_dict[DEPT._13C] and self._key == "d" and event.button == MouseButton.LEFT:
            x = event.xdata
            xmin, xmax = self._ax_dict[DEPT._13C].get_xlim()
            xscale: float = abs(xmax-xmin)*0.02
            for key, ppm in self._Result[DEPT._13C].items():
                if x-xscale < ppm < x+xscale:
                    self._nHydrogens_Atom[key][1] = 1
                    if self._nHydrogens_Atom[key][0] == 3:
                        self._nHydrogens_Atom[key][0] = 0
                    else:
                        self._nHydrogens_Atom[key][0] += 1

            self.clear_local_axes()
            self.draw_curve()
            self.draw_title()
            self.draw_x_axis()
            self.draw_status()
            self.draw_carbon_number()
            self._fig.canvas.draw_idle()

    def on_mouse_motion(self, event: MouseEvent) -> None:
        from matplotlib.backend_bases import MouseButton

        if event.inaxes is not None and event.inaxes is self._ax_dict.values() and event.button is MouseButton.LEFT:
            toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore

            if self._key == "t" and toolbar_mode == "zoom rect":
                self.draw_title()
                self._key = ""
                for x in self._ax_dict.values():
                    x.set_navigate_mode("ZOOM")
                self._fig.canvas.draw_idle()

            if self._key == "t" and toolbar_mode != "zoom rect" and event.ydata is not None:
                if event.inaxes == self._ax_dict[DEPT._90]:
                    self._thr[DEPT._90] = event.ydata
                elif event.inaxes == self._ax_dict[DEPT._135] and event.ydata >= 0:
                    self._thr[DEPT._135] = event.ydata
                elif event.inaxes == self._ax_dict[DEPT._135_down] and event.ydata < 0:
                    self._thr[DEPT._135_down] = event.ydata
                elif event.inaxes == self._ax_dict[DEPT._13C]:
                    self._thr[DEPT._13C] = event.ydata

                self.clear_local_axes()
                self.draw_curve()
                self.draw_title()
                self.draw_x_axis()
                self.draw_status()

                if event.inaxes in self._ax_dict.values() and self._key == "t" and event.inaxes == self._ax_selected:
                    start, end = event.inaxes.get_xlim()
                    event.inaxes.hlines(event.ydata, end, start, colors='k', linestyles="dashed", linewidth=1)  # type: ignore # nopep8

                self._fig.canvas.draw_idle()

    def draw_title(self) -> None:
        """Draw title on the plot."""
        self._ax_dict[DEPT._90].set_title(self._title, loc="left")

    def _show_help(self) -> None:
        pass
        """Display help information."""
        print("Help mode ")
        self._title = ("Press 'q' to Quit, 's' to Save file.\n"
                       "Press 't' to set the threshold.\n"
                       "Press 'd' to display the Carbon.\n")
        self.draw_title()
        self._fig.canvas.draw_idle()

    def connect(self) -> None:
        """Connect all event handlers."""
        self._cID_key: int = self._fig.canvas.mpl_connect(
            'key_press_event', self.on_key_press)  # type: ignore
        self._cID_button_press: int = self._fig.canvas.mpl_connect(
            'button_press_event', self.on_button_press)  # type: ignore
        self._cID_button_motion: int = self._fig.canvas.mpl_connect(
            'motion_notify_event', self.on_mouse_motion)  # type: ignore
        self._cID_button_release: int = self._fig.canvas.mpl_connect(
            'button_release_event', self.on_button_release)  # type: ignore

    def disconnect(self) -> None:
        self._fig.canvas.mpl_disconnect(self._cID_key)
        self._fig.canvas.mpl_disconnect(self._cID_button_press)
        self._fig.canvas.mpl_disconnect(self._cID_button_release)
        self._fig.canvas.mpl_disconnect(self._cID_button_motion)


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

    import platform
    _system: str = platform.system()
    if _system == "Linux":
        if args.dir is None:
            directory: Path = Path(
                "/home/vitamin/Simulation/38.Ergocalciferol(Vitamin_D2)/00.Spectra/bmse000510/nmr/set01")
        else:
            directory: Path = Path(args.dir)
    elif _system == "Darwin":
        if args.dir is None:
            directory: Path = Path(
                "/Users/chengwen-cheng/Desktop/Simulation/bmse000510/nmr/set01")
        else:
            directory: Path = Path(args.dir)
    else:
        print("  Only for ubuntu or Darwin system ...")
        print("  Exit and Close the program !!!")
        exit(0)

    # plot and indicate all peaks

    digrams: diagram = diagram(directory, args)
    digrams.connect()
    plt.ioff()
    plt.show()
    delete_all_files(pipe_fid_filename)


if __name__ == "__main__":
    main()
