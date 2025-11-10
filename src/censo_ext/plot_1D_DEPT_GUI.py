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
        "-i",
        "--input",
        dest="path",
        action="store",
        type=str,
        nargs="+",
        required=False,
        help="Provide the path of your pdata (under 1r folder) ",
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

    parser.add_argument(
        "--hidden",
        dest="hidden",
        action="store_true",
        help="Show the plot [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


# global variable
pipe_fid_filename = ".1d_pipe.fid"
peaks_fileName = "plot_1D_DEPT.peaks"


def Channel(args, path, thr: float, thr_ch3_180: float, ax: Axes, phase: float = 1.0) -> dict:
    dic, data = ng.bruker.read_pdata(str(path))
    udic = ng.bruker.guess_udic(dic, data)
    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(pipe_fid_filename)
    data = data.real*phase  # type: ignore
    uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)

    from censo_ext.Tools.spectra import numpy_thr_mean_3
    threshold: float = numpy_thr_mean_3(data.astype(np.float64))*thr
    if phase == 1:
        ax.hlines(threshold, args.end, args.start,
                  linestyles="dashdot", linewidth=0.5)
        if "DEPT_135" in str(path):
            threshold: float = numpy_thr_mean_3(
                data.astype(np.float64))*thr_ch3_180*(-1)
        ax.hlines(threshold, args.end, args.start,
                  linestyles="dashdot", linewidth=0.5)
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

    ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    ppm = np.linspace(ppm_1h_0, ppm_1h_1, data.shape[0])

    if isinstance(args.start, (int, float)) and isinstance(args.end, (int, float)):
        pass
    else:
        args.end, args.start = uc_1h.ppm_limits()

    from censo_ext.Tools.spectra import numpy_thr_mean_3
    threshold: float = 0
    if isinstance(thr, float):
        threshold: float = numpy_thr_mean_3(data.astype(np.float64))*thr
    # detect all peaks with a threshold
    from scipy.signal import find_peaks
    y_heighest = max(data)
    threshold += y_heighest * 0.01
    peaks, _ = find_peaks(data, height=threshold, width=1)

    # add markers for peak positions
    Result: dict[int, float] = dict()

    for n, peak in enumerate(peaks):
        ppm: float = uc_1h.ppm(peak)
        Result[n+1] = ppm

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

    if not args.hidden:
        print_arguments()

    path: dict[Path, Path] = {}

    ch1: Path = Path('13C')
    ch2: Path = Path('DEPT_90')
    ch3: Path = Path('DEPT_135')
    thr: dict[Path, float] = {ch1: 2.0, ch2: 20.0, ch3: 2.0}
    thr_ch3_180: float = 2.0

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
    fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
    ax: list[Axes] = fig.subplots(3, 1, sharex=True)  # type: ignore

    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                        top=0.90, wspace=0.05, hspace=0.05)

    channel: Path = ch1
    Result_ch1 = Channel(
        args, path=path[channel], thr=thr[channel], thr_ch3_180=thr_ch3_180, ax=ax[2])

    StAtoms: dict[int, int] = {key: -1 for key, value in Result_ch1.items()}

    channel: Path = ch2
    Result_ch2 = Channel(
        args, path=path[channel], thr=thr[channel], thr_ch3_180=thr_ch3_180, ax=ax[0])

    print("DEPT90             ppm")
    for idx1, ppm in Result_ch2.items():
        print(f"{idx1:6d} {ppm:>15.5f}")

    channel: Path = ch3
    Result_ch3 = Channel(
        args, path=path[channel], thr=thr[channel], thr_ch3_180=thr_ch3_180, ax=ax[1])

    Compare_two_dict(Result_ch1, Result_ch3, StAtoms, Label=3)
    Compare_two_dict(Result_ch1, Result_ch2, StAtoms, Label=1)
    print("DEPT135(up)        ppm")
    for idx1, ppm in Result_ch3.items():
        print(f"{idx1:6d} {ppm:>15.5f}")

    channel: Path = ch3
    Result_ch3_180 = Channel(args, path=path[channel], thr=thr[channel], thr_ch3_180=thr_ch3_180, ax=ax[2], phase=-1.0)  # nopep8
    Compare_two_dict(Result_ch1, Result_ch3_180, StAtoms, Label=2)

    print("DEPT135(down)      ppm")
    #
    for idx1, ppm in Result_ch3_180.items():
        print(f"{idx1:6d} {ppm:>15.5f}")
    # ic(StAtoms)

    for key, value in StAtoms.items():
        if value == -1:
            StAtoms[key] = 0

    channel: Path = ch1
    dic, data = ng.bruker.read_pdata(str(path[channel]))
    udic = ng.bruker.guess_udic(dic, data)
    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    pipe_fid_filename = ".1d_pipe.fid"
    ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(pipe_fid_filename)
    data = data.real  # type: ignore
    uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)

    # end ---------+--------- start
    # args.end                args.start
    # ppm_1h_0                ppm_1h_1

    ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    ppm: npt.NDArray[np.float64] = np.linspace(
        ppm_1h_0, ppm_1h_1, data.shape[0])

    if isinstance(args.start, (int, float)) and isinstance(args.end, (int, float)):
        pass
    else:
        args.end, args.start = uc_1h.ppm_limits()

    threshold: float = numpy_thr_mean_3(
        data.astype(np.float64))*thr[channel]

    # detect all peaks with a threshold
    from scipy.signal import find_peaks
    y_heighest = max(data)
    # y_lowest = min(data)
    threshold += float(y_heighest) * 0.01
    peaks, _ = find_peaks(data, height=threshold, width=1)

    # print the final data
    print("#   ID             ppm    nHydrogens")
    for n, peak in enumerate(peaks):
        height = data[int(peak)]
        ppm = uc_1h.ppm(peak)
        print(f"{n+1:6d} {ppm:>15.5f}        {StAtoms[n+1]:>3d}")

    # save to file
    if args.save:
        with open(peaks_fileName, "w") as f:
            sys.stdout = f
            print("#   ID             ppm    nHydrogens")
            for n, peak in enumerate(peaks):
                height = data[int(peak)]
                ppm = uc_1h.ppm(peak)
                print(f"{n+1:6d} {ppm:>15.5f}        {StAtoms[n+1]:>3d}")
            sys.stdout = sys.__stdout__

    # add markers for peak positions
    for n, peak in enumerate(peaks):
        height = data[int(peak)]
        ppm = uc_1h.ppm(peak)
        if ppm < args.end and ppm > args.start:
            if StAtoms[n+1] == 0:
                ax[2].text(ppm.tolist(), height*1.20,
                           str("C"), ha="center", va="center", rotation=90)
            elif StAtoms[n+1] == 1:
                ax[2].text(ppm.tolist(), height*1.20,
                           "CH", ha="center", va="center", rotation=90)
            else:
                ax[2].text(ppm.tolist(), height*1.20,
                           "CH"+rf'$_{str(StAtoms[n+1])}$', ha="center", va="center", rotation=90)

    fig.suptitle(args.path, fontsize=12, y=0.98)
    fig.text(0.5, 0.04, "$\\delta$ / ppm",
             ha="center", fontsize=12)
    ax[2].tick_params(axis="x", which="both", bottom=True,
                      top=False, labelbottom=True, labelsize=12)
    ax[2].spines["bottom"].set_visible(True)

    plt.show()
    delete_all_files(pipe_fid_filename)


if __name__ == "__main__":
    main()
