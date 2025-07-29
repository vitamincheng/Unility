#! /usr/bin/env python3
from ast import arguments
from pathlib import Path
from icecream import ic
import nmrglue as ng
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from sys import argv as sysargv
import argparse
import os
import sys

# global variable
peaks_fileName = "plot_1D_DEPT.peaks"

descr = """
________________________________________________________________________________
|                                          [10.05.2024] vitamin.cheng@gmail.com
| For Plot 1D sepctra in experiments using nmrglue module
| Usage: plot_1D_DEPT.py <geometry> [options]
| [Options]
| Input    : -i the pdata path(under 1r folder) [required]
|          : -start start point of chemical shift [default from data]
|          : -end   end point of chemical shift [default from data]
| Save     : --save saved the report of carbon [default false]
| Hidden   : -h show the plot [default False]
| Package  : Tools 
| Module   : spectra.py
|______________________________________________________________________________
"""
useit = """
    End     Endremove    Startremove                 Start
    +               +    +                               +
    +---------------+----+-------------------------------+
    lower field                               higher field
                        delta /ppm
    """


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        #        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
        add_help=False
    )  # argparse.RawDescriptionHelpFormatter) #,

    parser.add_argument(
        "-i",
        "--input",
        dest="path",
        action="store",
        type=str,
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
        default=False,
        help="Saved the report of carbon [default False]",
    )

    parser.add_argument(
        "-h",
        "--hidden",
        dest="hidden",
        action="store_true",
        default=False,
        help="Show the plot [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def Channel(args, path: dict, channel: Path, thr: float, phase: float = 1.0) -> dict[int, float]:
    dic, data = ng.bruker.read_pdata(str(path[channel]))
    udic = ng.bruker.guess_udic(dic, data)
    # ic(phase)
    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    pipe_fid_filename = ".1d_pipe.fid"
    ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(pipe_fid_filename)
    data = data.real*phase  # type: ignore
    uc_1h = ng.pipe.make_uc(dic, data)

    # end ---------+--------- start
    # args.end                args.start
    # ppm_1h_0                ppm_1h_1

    ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    ppm = np.linspace(ppm_1h_0, ppm_1h_1, data.shape[0])
    if args.start == None or args.end == None:
        args.end, args.start = uc_1h.ppm_limits()

    output = np.vstack((ppm, np.real(data))).T[::-1]
    # np.savetxt("output.dat", output, fmt=" %12.5f  %12.5e")
    from censo_ext.Tools.spectra import numpy_thr_mean_3
    threshold: float = 0
    if isinstance(thr, float):
        threshold: float = numpy_thr_mean_3(data.astype(np.float64))*thr
        # ic(threshold, thr)
    # detect all peaks with a threshold
    from scipy.signal import find_peaks
    y_heighest = max(data)
    y_lowest = min(data)
    # ic(len(data))
    # for H 65536 for C 131072
    # DEPT 90 32768 DEPT 32768
    threshold += y_heighest * 0.01
    peaks, _ = find_peaks(data, height=threshold, width=1)

    # add markers for peak positions
    Result: dict = dict()

    for n, peak in enumerate(peaks):
        ppm: float = uc_1h.ppm(peak)
        Result[n+1] = ppm
    return Result


def Compare_two_dict(CH1: dict, CH2: dict, StAtoms: dict, Label: int) -> None:
    from censo_ext.Tools.spectra import find_nearest
    for x in CH2.values():
        nearest_peak, idx0 = find_nearest(list(CH1.values()), x)
        # ic(nearest_peak-x)
        if (nearest_peak-x) < 0.05:
            # if StAtoms[idx0+1] == -1:
            StAtoms[idx0+1] = Label
            # else:
            #    print("some peaks is overlap")
            #    ic(StAtoms)
            #    ic(idx0+1)
            #    exit(0)
        else:
            print("some peaks is more than 0.02 ppm")
            ic()
            exit(0)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    # read in the Bruker data
    if args == argparse.Namespace():
        args = cml(descr)
    if not args.hidden:
        print(descr)  # Program description
        print("    provided arguments: {}".format(" ".join(sysargv)))

    path: dict[Path, Path] = {}

    ch1: Path = Path('13C')
    ch2: Path = Path('DEPT_90')
    ch3: Path = Path('DEPT_135')
    thr: dict[Path, float] = {ch1: 2.0, ch2: 20.0, ch3: 2.0}
    thr_ch3_180: float = 2.0
    directory: Path = Path("../../bmse000510/nmr/set01")
    path[ch1] = directory / ch1 / Path("pdata/1")
    path[ch2] = directory / ch2 / Path("pdata/1")
    path[ch3] = directory / ch3 / Path("pdata/1")

    channel: Path = ch1
    Result_ch1: dict = Channel(args, path, channel, thr[channel])

    # print("13C  #             ppm")
    # for idx1, ppm in Result_ch1.items():
    # print(f"{idx1:6d}", f"{ppm:>15.5f}")
    StAtoms: dict[int, int] = {key: -1 for key, value in Result_ch1.items()}

    channel: Path = ch2
    Result_ch2: dict = Channel(args, path, channel, thr[channel])

    # from Tools.spectra import find_nearest
    # nearest_peak = find_nearest(list(Result_ch1.values()), Result_ch2[1])
    # ic(nearest_peak[0]-Result_ch2[1])

    print("DEPT90             ppm")
    for idx1, ppm in Result_ch2.items():
        print(f"{idx1:6d} {ppm:>15.5f}")

    channel: Path = ch3
    Result_ch3: dict = Channel(args, path, channel, thr[channel])
    Compare_two_dict(Result_ch1, Result_ch3, StAtoms, Label=3)
    Compare_two_dict(Result_ch1, Result_ch2, StAtoms, Label=1)
    print("DEPT135(up)        ppm")
    for idx1, ppm in Result_ch3.items():
        print(f"{idx1:6d} {ppm:>15.5f}")

    channel: Path = ch3
    Result_ch3_180: dict = Channel(args, path, channel, thr_ch3_180, phase=-1.0)  # nopep8
    Compare_two_dict(Result_ch1, Result_ch3_180, StAtoms, Label=2)

    print("DEPT135(down)      ppm")
    #
    for idx1, ppm in Result_ch3_180.items():
        print(f"{idx1:6d} {ppm:>15.5f}")
    # ic(StAtoms)

    for key, value in StAtoms.items():
        if value == -1:
            StAtoms[key] = 0

    # ic(StAtoms)
    # channel:
    #
    #
    channel = ch1
    dic: dict
    data: npt.NDArray
    dic, data = ng.bruker.read_pdata(str(path[channel]))

    udic: dict = ng.bruker.guess_udic(dic, data)

    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    ng.pipe.write(".1d_pipe.fid", *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(".1d_pipe.fid")
    data = data.real  # type: ignore
    uc_1h = ng.pipe.make_uc(dic, data)

    # end ---------+--------- start
    # args.end                args.start
    # ppm_1h_0                ppm_1h_1

    ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    ppm: npt.NDArray[np.float64] = np.linspace(
        ppm_1h_0, ppm_1h_1, data.shape[0])
    if args.start == None or args.end == None:
        args.end, args.start = uc_1h.ppm_limits()

    output: npt.NDArray[np.float64] = np.vstack((ppm, np.real(data))).T[::-1]
    np.savetxt("output.dat", output, fmt=" %12.5f  %12.5e")
    from censo_ext.Tools.spectra import numpy_thr_mean_3
    threshold: float = numpy_thr_mean_3(
        data.astype(np.float64))*thr[channel]

    # detect all peaks with a threshold
    from scipy.signal import find_peaks
    y_heighest = max(data)
    y_lowest = min(data)
    # ic(len(data))
    # for H 65536 for C 131072
    # DEPT 90 32768 DEPT 32768
    threshold += y_heighest * 0.01
    peaks, _ = find_peaks(data, height=threshold, width=1)

    # plot and indicate all peaks
    fig = plt.figure(figsize=(11.7, 8.3), dpi=100)
    ax = fig.subplots()
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                        top=0.90, wspace=0.05, hspace=0.05)

    # print the final data
    print("#   ID             ppm    nHydrogens")
    for n, peak in enumerate(peaks):
        height = data[int(peak)]
        ppm = uc_1h.ppm(peak)
        print(f"{n+1:6d} {ppm:>15.5f}        {StAtoms[n+1]:>3d}")

    # save to file
    if args.save == True:
        original_stdout = sys.stdout
        with open(peaks_fileName, "w") as f:
            sys.stdout = f
            print("#   ID             ppm    nHydrogens")
            for n, peak in enumerate(peaks):
                height = data[int(peak)]
                ppm = uc_1h.ppm(peak)
                print(f"{n+1:6d} {ppm:>15.5f}        {StAtoms[n+1]:>3d}")
        sys.stdout = original_stdout

    # add markers for peak positions
    for n, peak in enumerate(peaks):
        height = data[int(peak)]
        ppm = uc_1h.ppm(peak)
        if ppm < args.end and ppm > args.start:
            # ax.scatter(ppm, height, marker="o", color="r",
            #           s=100, alpha=0.5)  # type: ignore
            # ax.text(ppm, height + threshold*5, f"{contour_heights[n]:12.3f}" , ha="center", va="center",rotation=90)
            if StAtoms[n+1] == 0:
                ax.text(ppm.tolist(), height*1.10,
                        str("C"), ha="center", va="center", rotation=90)
            elif StAtoms[n+1] == 1:
                ax.text(ppm.tolist(), height*1.10,
                        "CH", ha="center", va="center", rotation=90)
            else:
                ax.text(ppm.tolist(), height*1.10,
                        "CH"+r'$_{}$'.format(str(StAtoms[n+1])), ha="center", va="center", rotation=90)

            # ax.text(ppm.tolist(), height*1.05,
            #        str(n+1), ha="center", va="center")

    plt.plot(uc_1h.ppm_scale(), data, 'b', linewidth=1)
    plt.hlines(threshold, args.end, args.start,
               linestyles="dashdot", linewidth=0.5)
    plt.xlim(args.end, args.start)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", which="both", bottom=True,
                   top=False, labelbottom=True, labelsize=12)
    ax.tick_params(axis="y", which="both", left=False,
                   right=False, labelleft=False)
    ax.get_yaxis().set_visible(False)
    # plt.xlim(ppm_1h_0,ppm_1h_1)

    if y_lowest*(-1) < y_heighest*0.2:
        plt.ylim(-0.05*y_heighest, 1.10*y_heighest)
    else:
        plt.ylim(1.10*y_lowest, 1.10*y_heighest)
        # plt.ylim(-0.01*y_heighest, 1.10*y_heighest)
    fig.suptitle(args.path, fontsize=12, y=0.98)
    fig.text(0.5, 0.04, "$\\delta$ / ppm", ha="center", fontsize=12)
    plt.show()


if __name__ == "__main__":
    main()
