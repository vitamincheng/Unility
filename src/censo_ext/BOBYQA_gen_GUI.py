#!/usr/bin/env python
import argparse
from matplotlib.backend_bases import KeyEvent, MouseEvent
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from censo_ext import BOBYQA_gen
from censo_ext.Tools.anmrfile import AD_Normal
from censo_ext.Tools.datfile import CensoDat, Peaks_npz, unit_conversion
from censo_ext.Tools.utility import IsExist_bool, print_arguments

descr = """
_______________________________________________________________________________
| For generate the orcaS-BOBYQA.out
| Usages    : BOBYQA_gen_GUI.py <geometry> [options]
| [options]
| File      : -i input dat/npz file [default 1r.npz and output.npz]
| Dir       : -d the location directory [default .]
| Delete    : -d --delete Delete specific cID peaks [default None]
| Start     : -s --start Start ppm of Chemical Shift [default from Data]
| End       : -e --end End ppm of Chemical Shift [default from Data]
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
        nargs=2,
        default=["1r.npz", "output.npz"],
        help="Provide input one npz file [default ['1r.npz','output.npz']]",
    )

    parser.add_argument(
        "-d",
        "--dir",
        dest="dir",
        action="store",
        type=str,
        default=".",
        help="Location Directory [default .]",
    )

    parser.add_argument(
        "--delete",
        dest="delete",
        action="store",
        type=int,
        nargs="+",
        default=None,
        help="Delete specific cID peaks [default None]",
    )

    parser.add_argument(
        "-s",
        "-start",
        dest="start",
        action="store",
        type=float,
        default=None,
        help="Start ppm of Chemical Shift [default from Data]",
    )

    parser.add_argument(
        "-e",
        "-end",
        dest="end",
        action="store",
        type=float,
        default=None,
        help="End ppm of Chemical Shift [default from Data]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


class diagram:
    """A class for creating and managing a spectrum diagram with interactive editing capabilities."""

    def __init__(self, args: argparse.Namespace, intensit: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]],
                 uc: tuple[unit_conversion, unit_conversion], peaks_npz) -> None:
        """Initialize the diagram with spectrum data and plotting setup.

        Args:
            fileName: Path to the spectrum file
            peaks_npz: Peaks data container
            intensit: Intensity array of the spectrum
            uc: Unit conversion object for ppm scale
            thres: Threshold value for spectrum
            ng_1r_peaks: Array of peaks information
        """
        self._fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
        self._ax_list: list[Axes] = self._fig.subplots(
            2, 1, sharex=True)  # type: ignore
        self._ax_top: Axes = self._ax_list[0]
        self._ax_bottom: Axes = self._ax_list[1]
        self._fig.subplots_adjust(left=0.07, right=0.93, bottom=0.10,
                                  top=0.90, wspace=0.05, hspace=0.05)

        self._1r_fileName: Path = Path(args.file[0])
        self._output_fileName: Path = Path(args.file[1])
        self._peaks_npz: Peaks_npz = peaks_npz
        self._intensit_top: npt.NDArray[np.float64] = intensit[0]
        self._intensit_bottom: npt.NDArray[np.float64] = intensit[1]
        self._button_xy: tuple[float, float] | None = None
        self._uc_top: unit_conversion = uc[0]
        self._uc_bottom: unit_conversion = uc[1]
        self._top_status_int: list[int] = peaks_npz.get_cIDs_center_peaks()[0]
        self._bottom_status_select_int: list[int] = []
        self._bottom_status_delete_int: list[int] = []
        self._bottom_status_str: str
        self._title_bottom: str = "Edit mode.\nPress 'h' to help. "
        self.draw_title()
        self.draw_status()
        self._key: str = ""
        self._ax_selected = None
        self._bottom_selected: bool = False
        self._bottom_selected_limits = None | tuple[float, float]
        self._AD_orcaS: AD_Normal = AD_Normal()
        self.method_load_orcaS()
        self.x_boundary(args, uc)

    def x_boundary(self, args: argparse.Namespace, uc: tuple[unit_conversion, unit_conversion]) -> None:
        if args.start is not None and args.end is not None and args.start < args.end:
            self._start: float = args.start
            self._end: float = args.end
        else:
            self._start = uc[1].ppm_limits()[0] if uc[1].ppm_limits(
            )[0] > uc[0].ppm_limits()[0] else uc[0].ppm_limits()[0]
            self._end = uc[1].ppm_limits()[1] if uc[1].ppm_limits(
            )[1] < uc[0].ppm_limits()[1] else uc[0].ppm_limits()[1]

    def method_load_orcaS(self) -> None:
        if (self._AD_orcaS.Exist()):
            self._AD_orcaS.method_load_files()
            if isinstance(self._AD_orcaS.ChemicalShifts, dict):
                self.OrcaS: npt.NDArray[np.float64] = np.array(
                    list(self._AD_orcaS.ChemicalShifts.items()))
                print("\n  ===== Loading data OrcaS.out of Average Directory =====")
                print(f"{self._AD_orcaS._orcaS}\n{self.OrcaS}")
            else:
                print("  The OrcaS.out in Averaage Directory is not dict format !!!")
                print("  Exit and Close the program !!!")
                exit(0)

    def on_key_press(self, event: KeyEvent) -> None:
        """Callback function for key press events."""
        if event.key == 'q':
            print("Quitting the application.")
            plt.close(event.canvas.figure)
        elif event.key == 'enter':
            self._execute_action()

        elif event.key == 'escape':
            self._reset_to_edit_mode()

        elif event.key == 'h':
            self._show_help()

        elif event.key == 'g':
            self._generate_file()

        elif event.key == 'd':
            self._set_delete_mode()

        elif event.key == 'e':
            self._set_select_mode()

        elif event.key == 'f':
            self._redraw_full()

    def _execute_action(self) -> None:
        """Execute the current action based on key mode."""

        self._bottom_status_select_int = list(
            map(int, set(self._bottom_status_select_int)))
        self._bottom_status_select_int.sort()

        if self._key == "d":
            for x in self._bottom_status_delete_int:
                if x in self._bottom_status_select_int:
                    self._bottom_status_select_int.remove(x)
            self._bottom_status_delete_int = []

        self.clear_local_axes()
        self.draw_curve()
        self.draw_title()
        self.draw_status()
        self.draw_scatter_numbers()
        self._fig.canvas.draw_idle()

    def clear_local_axes(self) -> None:

        top_x, top_y = self._ax_top.get_xlim(), self._ax_top.get_ylim()
        bottom_x, bottom_y = self._ax_bottom.get_xlim(), self._ax_bottom.get_ylim()
        self._ax_top.clear()
        self._ax_bottom.clear()
        self._ax_top.set_xlim(top_x)
        self._ax_top.set_ylim(top_y)
        self._ax_bottom.set_xlim(bottom_x)
        self._ax_bottom.set_ylim(bottom_y)

    def _reset_to_edit_mode(self) -> None:
        """Reset to edit mode."""
        print("Edit mode : ")
        self._key = 'escape'

        self.clear_local_axes()
        self._title_bottom = "Edit mode.\nPress 'h' to help. "
        self.draw_curve()
        self.draw_title()
        self.draw_status()
        self.draw_scatter_numbers()
        self._fig.canvas.draw_idle()

    def _show_help(self) -> None:
        """Display help information."""
        print("Help mode : ")
        self._bottom_status_select_int = []
        self._title_bottom = ("Press 'q' to Quit, 'g' to Generate file.\n"
                              "Press 'd' to Delete, 'f' to Full screen\n"
                              "'e' sElect mode, 'Esc' retrun to Edit mode\n"
                              "'Enter' execute mode, ")
        self.draw_title()
        self._fig.canvas.draw_idle()

    def _generate_file(self) -> None:
        """Save the current peaks data."""
        self._key = 'g'
        print("Generate mode : ", end="")
        self._title_bottom = "Generate file : "

        args_x: dict = {"file": "peaks.npz", "comb": 100, "start": None, "end": None,
                        "index": self._bottom_status_select_int, "delete": None}

        BOBYQA_gen.main(argparse.Namespace(**args_x))

        self.draw_title()
        self._fig.canvas.draw_idle()

    def _set_delete_mode(self) -> None:
        """Set edit mode based on key pressed."""
        self._key = 'd'
        print("Delete mode : ", end="")
        self._title_bottom = "Delete mode : "
        self.draw_title()
        self._fig.canvas.draw_idle()

    def _set_select_mode(self) -> None:
        """Set select mode."""
        toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore
        if toolbar_mode == "zoom rect":
            self._ax_top.set_navigate_mode(None)
        else:
            self._key = 'e'
            print("Select mode : ", end="")
            self._title_bottom = "Select mode : "
            self.draw_title()
            self._fig.canvas.draw_idle()

    def _redraw_full(self) -> None:
        """Redraw the full spectrum."""

        xmin, xmax, ymin, ymax = plt.axis()
        self._ax_top.clear()
        self._ax_bottom.clear()
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

        self.draw_x_axis()
        self.draw_curve()
        self.draw_title()
        self.draw_status()
        self.draw_scatter_numbers()
        self._fig.canvas.draw_idle()

    def on_button_release(self, event: MouseEvent) -> None:
        """Handle button release events."""
        if event.inaxes is self._ax_selected and self._button_xy is not None:

            if event.xdata is None or event.ydata is None:
                release_x: float = self._button_xy[0]
                release_y: float = self._button_xy[1]
            else:
                release_x: float = float(event.xdata)
                release_y: float = float(event.ydata)

            distance: np.float64 = np.sqrt((release_x-self._button_xy[0]) **
                                           2 + (release_y-self._button_xy[1])**2)
            tolerance = 0.01

            if self._key == 'e':
                if event.inaxes is self._ax_bottom:
                    self._bottom_selected = True
                    self._bottom_selected_limits = self._button_xy[0], release_x
                    end, start = self._bottom_selected_limits
                    a: list = []
                    for x in self.OrcaS:
                        if start < x[1] < end:
                            a.append(int(x[0]))
                    self._bottom_status_select_int = a

                    self.clear_local_axes()
                    self.draw_curve()
                    self.draw_title()
                    self.draw_status()
                    self.draw_scatter_numbers()
                    self._fig.canvas.draw_idle()

                    return None
            cID: int
            if distance < tolerance:
                if self._key == "d":
                    for x in self.OrcaS:
                        if release_x-0.020 < x[1] < release_x+0.020:
                            cID = int(x[0])
                            print(f" {cID}", end="")
                            self._bottom_status_delete_int.append(cID)
            else:
                return None

            self._button_xy = None
            self.draw_status()
            self._fig.canvas.draw_idle()

    def on_button_press(self, event: MouseEvent) -> None:
        """Handle button press events."""
        from matplotlib.backend_bases import MouseButton

        toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore
        if self._key == "e" and toolbar_mode == "zoom rect":
            self._ax_bottom.set_navigate_mode("ZOOM")
            self._ax_top.set_navigate_mode("ZOOM")

        elif (event.inaxes == self._ax_bottom) and event.button == MouseButton.LEFT:
            self._button_xy = event.xdata, event.ydata  # type: ignore
            self._ax_selected = event.inaxes

    def on_mouse_motion(self, event: MouseEvent) -> None:
        """Handle mouse motion events."""
        from matplotlib.backend_bases import MouseButton

        if event.inaxes is self._ax_selected and event.button is MouseButton.LEFT:
            toolbar_mode = self._fig.canvas.manager.toolbar.mode  # type: ignore

            if self._key == "e" and toolbar_mode == "zoom rect":
                self._title_bottom = ""
                self.draw_title()
                self._key = ""
                self._ax_bottom.set_navigate_mode("ZOOM")
                self._ax_top.set_navigate_mode("ZOOM")
                self._fig.canvas.draw_idle()
            if self._key == "e" and toolbar_mode != "zoom rect":
                bottom_y = self._ax_bottom.get_ylim()

                self.clear_local_axes()
                self.draw_curve()
                self.draw_title()
                self.draw_status()
                self.draw_scatter_numbers()

                if event.inaxes is self._ax_bottom:
                    self._ax_bottom.vlines(self._button_xy[0], self._button_xy[1]-bottom_y[0]*0.01, self._button_xy[1]+bottom_y[0]*0.01, colors='k', linestyles="solid", linewidth=2)  # type: ignore # nopep8
                    self._ax_bottom.vlines(event.xdata, self._button_xy[1]-bottom_y[0]*0.01, self._button_xy[1]+bottom_y[0]*0.01, colors='k', linestyles="solid", linewidth=2)  # type: ignore # nopep8
                    self._ax_bottom.hlines(self._button_xy[1], self._button_xy[0], event.xdata, colors='k', linestyles="solid", linewidth=1)  # type: ignore # nopep8

                self._fig.canvas.draw_idle()

    def connect(self) -> None:
        """Connect all event handlers."""
        self._cID_key = self._fig.canvas.mpl_connect(
            'key_press_event', self.on_key_press)  # type: ignore
        self._cID_button_press = self._fig.canvas.mpl_connect(
            'button_press_event', self.on_button_press)  # type: ignore
        self._cID_button_motion = self._fig.canvas.mpl_connect(
            'motion_notify_event', self.on_mouse_motion)  # type: ignore
        self._cID_button_release = self._fig.canvas.mpl_connect(
            'button_release_event', self.on_button_release)  # type: ignore

    def disconnect(self) -> None:
        """Disconnect all event handlers."""
        self._fig.canvas.mpl_disconnect(self._cID_key)
        self._fig.canvas.mpl_disconnect(self._cID_button_press)
        self._fig.canvas.mpl_disconnect(self._cID_button_release)
        self._fig.canvas.mpl_disconnect(self._cID_button_motion)

    def draw_status(self) -> None:
        """Draw status on the plot."""
        self._ax_top.set_title(
            f"Select : {self._top_status_int}\nSizes of Select : {len(self._top_status_int)}", loc="right", fontsize=10)
        self._ax_bottom.set_title(
            f"Select : {self._bottom_status_select_int}\nSizes of Select : {len(self._bottom_status_select_int)}\n\
              Delete : {self._bottom_status_delete_int}", loc="right", y=0, fontsize=10)
        # self._fig.canvas.draw_idle()

    def draw_title(self) -> None:
        """Draw title on the plot."""
        self._ax_bottom.set_title(self._title_bottom, loc="left", y=0)
        # self._fig.canvas.draw_idle()

    def draw_scatter_numbers(self) -> None:
        """Draw integral curves on the plot."""
        Data: list[tuple[int, npt.NDArray[np.float64], npt.NDArray[np.float64]]
                   ] = self._peaks_npz.method_integrate(self._intensit_top)
        for cID, peak_int, peak_scale in Data:
            self._ax_top.text(peak_scale[0], 0.5 * peak_int.sum() / 100./4 + peak_int.max()*0.8, str(cID),
                              fontsize=8)
        a: npt.NDArray[np.float64] = self._peaks_npz.get_cIDs_center_peaks()
        for ppm in a[1]:
            index: int = self._uc_top.index(ppm)
            height = float(self._intensit_top[index])
            self._ax_top.scatter(ppm, height, marker="o",
                                 color="r", s=30, alpha=0.5)

        if len(self._bottom_status_select_int) != 0:
            y_highest, _ = self._ax_bottom.get_ylim()
            for key, ppm in self._AD_orcaS.ChemicalShifts.items():  # type: ignore
                if key in self._bottom_status_select_int:
                    index: int = self._uc_bottom.index(ppm)
                    height = float(self._intensit_bottom[index])
                    self._ax_bottom.scatter(ppm, height, marker="o",
                                            color="y", s=30, alpha=0.5)
                    self._ax_bottom.text(ppm, height+y_highest*0.02, str(key),
                                         fontsize=8)

    def draw_curve(self) -> None:
        """Draw the main spectrum curve."""
        self._ax_top.plot(self._uc_top.ppm_scale(),
                          self._intensit_top, 'b', linewidth=1)
        self._ax_bottom.plot(self._uc_bottom.ppm_scale(),
                             self._intensit_bottom, 'k', linewidth=1)

    def draw_x_axis(self) -> None:
        """Draw x-axis configuration."""
        y_heighest_1r = float(np.max(self._intensit_top))
        y_lowest_1r = float(np.min(self._intensit_top))
        y_heighest_output = float(np.max(self._intensit_bottom))
        y_lowest_output = float(np.min(self._intensit_bottom))
        self._ax_bottom.set_xlim(self._end, self._start)

        # top
        self._ax_top.spines["right"].set_visible(False)
        self._ax_top.spines["bottom"].set_visible(False)
        self._ax_top.spines["top"].set_visible(False)
        self._ax_top.spines["left"].set_visible(False)
        self._ax_top.get_yaxis().set_visible(False)
        self._ax_top.get_xaxis().set_visible(False)

        # bottom
        self._ax_bottom.spines["right"].set_visible(False)
        self._ax_bottom.spines["bottom"].set_visible(True)
        self._ax_bottom.spines["left"].set_visible(False)
        self._ax_bottom.spines["top"].set_visible(False)
        self._ax_bottom.tick_params(axis="x", which="both", top=False,
                                    bottom=True, labelbottom=True, labelsize=12)
        # self._ax_output.tick_params(axis="y", which="both", left=False,
        #                            right=False, labelleft=False)
        self._ax_bottom.get_yaxis().set_visible(False)

        # If phase is -1, it will adjust the y axis
        if y_lowest_1r*(-1) < y_heighest_1r*0.2:
            self._ax_top.set_ylim((-0.05*y_heighest_1r, 1.10*y_heighest_1r))
        else:
            self._ax_top.set_ylim((1.10*y_lowest_1r, 1.10*y_heighest_1r))

        if y_lowest_output*(-1) < y_heighest_output*0.2:
            self._ax_bottom.set_ylim(
                (1.10*y_heighest_output, -0.05*y_heighest_output))
        else:
            self._ax_bottom.set_ylim(
                (1.10*y_heighest_output, 1.10*y_lowest_output))

        self._fig.suptitle(str(self._1r_fileName), fontsize=12, y=0.98)
        self._fig.text(0.5, 0.04, "$\\delta$ / ppm", ha="center", fontsize=12)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    plt.rcParams['keymap.save'].remove('s')
    plt.rcParams['keymap.fullscreen'].remove('f')
    plt.rcParams['keymap.back'].remove('c')
    plt.rcParams['keymap.grid'].remove('g')
    plt.rcParams['toolbar'] = 'toolbar2'
    plt.ion()

    args_file0: Path = Path(args.dir) / Path(args.file[0])
    args_file1: Path = Path(args.dir) / Path(args.file[1])

    if not IsExist_bool(args_file0) or not IsExist_bool(args_file1):
        return

    # Load the data from dat/npz file
    # 1r.npz of file[0]     Data_0     ppm[0]       top of diagram
    # output.npz of file[1] Data_1     ppm[1]    bottom of diagram
    censo_0: CensoDat = CensoDat(args_file0)
    censo_1: CensoDat = CensoDat(args_file1)

    Data_0 = censo_0.get_Dat().T
    Data_1 = censo_1.get_Dat().T

    ppm: tuple[npt.NDArray[np.float64],
               npt.NDArray[np.float64]] = Data_0[0], Data_1[0]
    intensit: tuple[npt.NDArray[np.float64],
                    npt.NDArray[np.float64]] = Data_0[1], Data_1[1]

    uc: tuple[unit_conversion, unit_conversion] = unit_conversion(
        ppm[0]), unit_conversion(ppm[1])
    peaks_npz: Peaks_npz = Peaks_npz(uc[0])

    peaks_npz.method_read_file()
    print("  ========== Before ==========")
    peaks_npz.method_print()

    diagrams: diagram = diagram(
        args, intensit, uc, peaks_npz)

    diagrams.draw_scatter_numbers()
    diagrams.draw_curve()
    diagrams.draw_x_axis()
    diagrams.connect()

    plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()
