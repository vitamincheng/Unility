#! /usr/bin/env python3
from icecream import ic
import matplotlib.pyplot as plt
from nmrglue.fileio.fileiobase import unit_conversion
from censo_ext.Tools.spectra import numpy_threshold_3, numpy_threshold_10
from sys import argv as sysargv
import numpy as np
import nmrglue as ng
import argparse
import os

descr = """
________________________________________________________________________________
|                                          [10.05.2024] vitamin.cheng@gmail.com
| For Plot 2D sepctra in experiments using nmrglue module
| Usage: plot_2D_exp.py <geometry> [options]                  
| [Options]
| Input    : -i the pdata path(under 2rr folder) [required] 
|          : -F1 F1 shift, usually is Carbon   [default 0.0] 
|          : -F2 F2 shift, usually is Hydrogen [default 0.0]
| Hidden   : -h show the plot [default False] 
| Package  : Tools 
| Module   : spectra.py
|______________________________________________________________________________
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
        help="Provide the path of your pdata (under 2rr folder) ",
    )

    parser.add_argument(
        "-F2",
        "--DeltaF2",
        dest="DeltaF2",
        action="store",
        type=float,
        required=False,
        default=0,
        help="Chemical Shift for shift in F2 [default 0.0]",
    )

    parser.add_argument(
        "-F1",
        "--DeltaF1",
        dest="DeltaF1",
        action="store",
        type=float,
        required=False,
        default=0.0,
        help="Chemical Shift for shift in F1 [default 0.0]",
    )

    parser.add_argument(
        "-h",
        "--hidden",
        dest="hidden",
        action="store",
        type=bool,
        required=False,
        default=False,
        help="Show the plot [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def read_from_bruker(fileName, DeltaF1, DeltaF2):

    import nmrglue as ng
    dic, data_bruker = ng.bruker.read_pdata(fileName)
    udic = ng.bruker.guess_udic(dic, data_bruker)

    # ic(udic)
    udic[0]["car"] = udic[0]["car"] - DeltaF1*udic[0]['obs']
    udic[1]["car"] = udic[1]["car"] - DeltaF2*udic[1]['obs']
    import numpy as np
    data_bruker = data_bruker.astype(np.complex128)
    C = ng.convert.converter()
    C.from_bruker(dic, data_bruker, udic)
    pipe_fid_filename = ".2d_pipe.fid"
    ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
    udic, data = ng.pipe.read(pipe_fid_filename)
    data = data.real

    # Important
    # This Bruker is Quadruple Detector, added below two parameter (udic)
    # other optional is :
    # in ng.pipe.make_uc function in pipe.py
    # Remark  171 lines: #size = size / 2
    ###
    udic['FDF1QUADFLAG'] = 1.0
    udic['FDF2QUADFLAG'] = 1.0

    return udic, data


def read_udic(udic, data) -> tuple[unit_conversion, unit_conversion]:
    uc_1h: unit_conversion = ng.pipe.make_uc(udic, data, dim=1)
    uc_13c: unit_conversion = ng.pipe.make_uc(udic, data, dim=0)
    # ppm_1h = uc_1h.ppm_scale()
    # ppm_13c = uc_13c.ppm_scale()
    return uc_1h, uc_13c


def plot_2D_Basic(udic, data, uc_1h, uc_13c):
    # create the figure
    import numpy as np
    ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    ppm_13c_0, ppm_13c_1 = uc_13c.ppm_limits()
    fig = plt.figure(figsize=(11.7, 8.3), dpi=100)

    gs = fig.add_gridspec(2, 2,  width_ratios=(1, 19), height_ratios=(1, 9),
                          left=0.03, right=0.97, bottom=0.03, top=0.97,
                          wspace=0.1, hspace=0.1)

    ax = fig.add_subplot(gs[1, 1])
    ax_histx = fig.add_subplot(gs[0, 1], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 0], sharey=ax)
    ax_histx.get_xaxis().set_visible(False)
    ax_histx.get_yaxis().set_visible(False)
    ax_histx.axis('off')
    ax_histy.get_xaxis().set_visible(False)
    ax_histy.get_yaxis().set_visible(False)
    ax_histy.axis('off')

    x_axis_data = np.sum(data, axis=0)
    y_axis_data = np.sum(data, axis=1)

    # height_x_axis_data =np.max(x_axis_data)
    # height_y_axis_data =np.max(y_axis_data)

    ax_histx.plot(uc_1h.ppm_scale(), x_axis_data, color='k', linewidth=1)
    ax_histy.plot(-y_axis_data, uc_13c.ppm_scale(), color='k', linewidth=1)

    contour_thr = numpy_threshold_10(data)
    maximum = np.max(data)

    import matplotlib
    # type: ignore  # contour map (colors to use for contours)
    cmap = matplotlib.cm.Blues_r  # type: ignore nopep8
    contour_start: float = contour_thr  # contour level start value
    contour_num = 20                                # number of contour levels
    # scaling factor between contour levels
    contour_factor = 1.20

    import numpy as np
    cl = contour_start * contour_factor ** np.arange(contour_num)

    # plot the contours
    ax.contour(data, cl, cmap=cmap, extent=(
        ppm_1h_0, ppm_1h_1, ppm_13c_0, ppm_13c_1))

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_ylabel(udic['FDF1LABEL'] + " (ppm)")
    ax.set_xlabel(udic['FDF2LABEL'] + " (ppm)")
    ax.yaxis.set_label_coords(-0.10, 0.10)
    ax.xaxis.set_label_coords(0.95, 1.10)
    # ax.set_title("Protein 2D HSQC Spectrum")
    # ax.set_title("Protein 2D HMBC Spectrum")
    ax.set_xlim(ppm_1h_0, ppm_1h_1)
    ax.set_ylim(ppm_13c_0, ppm_13c_1)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    return ax


def cal_contour_peak(data, contour_thr_factor: float = 2):

    from skimage.morphology import extrema
    contour_maxima_thr = numpy_threshold_3(
        data) * contour_thr_factor + np.max(data)*0.01
    h_maxima = extrema.h_maxima(data, contour_maxima_thr)
    max_peaks: list = []
    for idy, y in enumerate(h_maxima):
        for idx, x in enumerate(y):
            if x == 1:
                max_peaks.append([idx, idy])
    return max_peaks


def main(args=argparse.Namespace()) -> None:
    import sys
    from icecream import ic
    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    # Important
    # Remember to fix the pipe.py in ng.pipe.make_uc to remark   171 lines: size = size / 2
    ###
    if args.path == None:
        args.path = "../../bmse000510/nmr/set01/1H_13C_HSQC/pdata/1/"
        # args.path = "../../bmse000510/nmr/set01/1H_13C_HMBC/pdata/1/"
        # args.path = "../../bmse000405/nmr/set01/HH_TOCSY/pdata/1/"
    udic, data = read_from_bruker(args.path, args.DeltaF1, args.DeltaF2)
    uc_1h, uc_13c = read_udic(udic, data)

    x_axis_data = np.sum(data, axis=0)
    y_axis_data = np.sum(data, axis=1)
    # ic(y_axis_data)
    from censo_ext.Tools.spectra import numpy_threshold_3
    y_thr: float = numpy_threshold_3(y_axis_data)
    x_thr: float = numpy_threshold_3(x_axis_data)
    # ic(y_thr)
    import matplotlib.axes as axes
    ax: axes.Axes | None = None
    if args.hidden == False:
        ax = plot_2D_Basic(udic, data, uc_1h, uc_13c)
    max_peaks: list = cal_contour_peak(data, contour_thr_factor=1)
    max_peaks = [a for a in list(
        max_peaks) if x_axis_data[a[0]] > x_thr and y_axis_data[a[1]] > y_thr]
    y_peaks: list = sorted(set([a[1] for a in list(max_peaks)]))
    # y_peaks: list = sorted(list(set(np.array(max_peaks).T[1])))
    # ic(y_peaks)
    x_grobal_maximum = x_axis_data.max()

    for y_idx in y_peaks:
        xslice = data[y_idx, :]
        maximum = xslice.max()
        x_sum = xslice.sum()
        xright = uc_1h.ppm(xslice.size)
        if args.hidden == False and ax != None:
            # ax.plot(uc_1h.ppm_scale(), -xslice/maximum *
            #        5 + uc_13c.ppm(y_idx), linewidth=0.5)
            ax.plot(uc_1h.ppm_scale(), -xslice/x_grobal_maximum *
                    5*4 + uc_13c.ppm(y_idx), linewidth=0.5)
            ax.text(xright, uc_13c.ppm(y_idx), f"{uc_13c.ppm(
                y_idx):12.3f}", ha="right", va="center", fontsize=6)
            for x in [x1 for x1, y1 in max_peaks if y1 == y_idx]:
                if data[y_idx][x] >= maximum * 0.5:
                    ax.scatter(uc_1h.ppm(x), uc_13c.ppm(y_idx), marker="o", color="r", s=300, alpha=0.5)  # type: ignore # nopep8
        # if y_idx == 52:
        #    ax.plot(uc_1h.ppm_scale(), -xslice/maximum *
        #            5 + uc_13c.ppm(y_idx), linewidth=0.5)
        #    ax.text(xright, uc_13c.ppm(y_idx), f"{uc_13c.ppm(
        #        y_idx):12.3f}", ha="right", va="center", fontsize=6)

    print("#          13C             1H       ")
    for x, y in max_peaks:
        xslice = data[y, :]
        maximum = xslice.max()
        if data[y][x] >= maximum * 0.5:
            print(f"{uc_13c.ppm(y):>15.4f} {uc_1h.ppm(x):>15.4f}")
            # ax.scatter(uc_1h.ppm(peak[0]), uc_13c.ppm(peak[1]), marker="o", color="r", s=300, alpha=0.5)  # type: ignore # nopep8

    from censo_ext.Tools.utility import save_figure
    save_figure()
    if args.hidden == False:
        plt.show()


if __name__ == "__main__":
    main()
