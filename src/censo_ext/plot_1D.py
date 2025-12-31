#! /usr/bin/env python
from nmrglue.fileio.fileiobase import unit_conversion
from censo_ext.Tools.utility import print_arguments
from icecream import ic
import nmrglue as ng
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import numpy as np
import numpy.typing as npt
import argparse

# global variable
peaks_fileName = "plot_1D.peaks"

descr = """
________________________________________________________________________________
| Usage: plot_1D.py <geometry> [options]
| For Plot 1D sepctra in experiments using nmrglue module
| For bruker file or jcamp-dx file and convert to dat/npz file
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
        dest="file",
        action="store",
        type=str,
        required=False,
        default=None,
        help="Provide the file of your jcamp-dx ",
    )

    parser.add_argument(
        "-d",
        "--dir",
        dest="path",
        action="store",
        type=str,
        required=False,
        default=None,
        help="Provide the path of your pdata (under 1r folder) ",
    )

    parser.add_argument(
        "-o",
        "--out",
        dest="out",
        action="store",
        type=str,
        required=False,
        help="Provide the filename of output files(dat/npz) ",
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
        "-thr",
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
        "--auto",
        dest="auto",
        action="store_true",
        help="Automatically integrate the peaks [default False]",
    )

    parser.add_argument(
        "--save",
        dest="save",
        action="store_true",
        help="Saved the automatically Integral of the peaks to plot_1D_peaks.out [default False]",
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

    parser.add_argument(
        "--show",
        dest="show",
        action="store_true",
        help="Show the spectra in screen [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    if args.path is None and args.file is None:

        # For test
        args.path = "../../Simulation/bmse000510/nmr/set01/1H/pdata/1"
        # args.path = "../../Simulation/bmse000510/nmr/set01/13C/pdata/1"
        # args.path = "../../Simulation/bmse000510/nmr/set01/DEPT_90/pdata/1"
        # args.path = "../../Simulation/bmse000510/nmr/set01/DEPT_135/pdata/1"

        # args.file = "../../Simulation/@@Exapmles/38.Ergocalciferol(Vitamin D2)/HMDB0000900_nmroned_1593_from_nmrml.jdx"

        # This is Normal Mode
        # print("  Only for one parameters in ars.path and args.file")
        # print("  Exit and Close the program !!!")
        # exit(0)

    if args.path is not None:
        try:
            dic, data = ng.bruker.read_pdata(args.path)
            udic: dict = ng.bruker.guess_udic(dic, data)
            C = ng.convert.converter()
            C.from_bruker(dic, data, udic)
            pipe_fid_filename: str = ".1d_pipe.fid"
            ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
            dic, data = ng.pipe.read(pipe_fid_filename)
        except OSError:
            try:
                dic, data = ng.varian.read(args.path)
                # ic(dic)
                tn = dic['procpar']['tn']['values']
                sw = float(dic['procpar']['sw']['values'][0])
                sfrq = float(dic['procpar']['sfrq']['values'][0])
                rp = float(dic['procpar']['rp']['values'][0])
                # ic(dic['procpar'])
                ic(tn, sw, sfrq, rp)
                udic: dict = ng.varian.guess_udic(dic, data)
                C = ng.convert.converter()
                C.from_varian(dic, data, udic)
                pipe_fid_filename: str = ".1d_pipe.fid"
                ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
                dic, data = ng.pipe.read(pipe_fid_filename)

                match tn:
                    case ['H1']:
                        dic['FDF2SW'] = dic['FDF2SW']*sw/sfrq    # for hydrogen
                    case ['C13']:
                        dic['FDF2SW'] = dic['FDF2SW']*sw/sfrq    # for carbon
                    case _:
                        print('  Only for H1 and C13')
                        exit(0)
                dic['FDF2ORIG'] = sfrq*rp
                ic(dic)
                # process the spectrum
                data = ng.proc_base.zf_size(data, 32768)
                data = ng.proc_base.fft(data)               # Fourier transform
                data = ng.proc_autophase.autops(
                    data=data, fn="acme")  # type: ignore
                # discard the imaginaries
                data = ng.proc_base.di(data)
                data = ng.proc_base.rev(data)               # reverse the data

            except OSError:
                print("  Only for bruker and varian ")
                exit(0)

    if args.file is not None:
        dic, data = ng.jcampdx.read(args.file)  # type: ignore
        dic = ng.jcampdx.guess_udic(dic, data)
        C = ng.convert.converter()
        C.from_universal(dic, data)
        pipe_fid_filename: str = ".1d_pipe.fid"
        ng.pipe.write(pipe_fid_filename, *C.to_pipe(), overwrite=True)
        dic, data = ng.pipe.read(pipe_fid_filename)

    data: npt.NDArray = data.real*args.phase  # type: ignore
    uc: unit_conversion = ng.pipe.make_uc(dic, data)  # type: ignore

    # end ---------+--------- start
    # args.end                args.start
    # limit_ppm[0]            limit_ppm[1]

    limit_ppm = uc.ppm_limits()
    ppm: npt.NDArray[np.float64] = np.linspace(
        limit_ppm[0], limit_ppm[1], data.shape[0])
    if args.start is None or args.end is None:
        args.end, args.start = uc.ppm_limits()

    if args.out:
        output: npt.NDArray[np.float64] = np.vstack(
            (ppm, np.real(data))).T[::-1]
        from censo_ext.Tools.utility import save_simulation_spectra_file
        save_simulation_spectra_file(args.out, output)

    from censo_ext.Tools.spectra import numpy_thr_mean_3
    threshold: float = numpy_thr_mean_3(data)*args.thr

    # detect all peaks with a threshold
    # from scipy.signal import find_peaks
    y_heighest: float = float(np.max(data))
    y_lowest: float = float(np.min(data))
    # ic(len(data))
    # for H 65536 for C 131072
    # DEPT 90 32768 DEPT 32768
    threshold += y_heighest * 0.01
    new_peaks = ng.peakpick.pick(
        data=data, pthres=threshold, algorithm="downward")
    ic(new_peaks)

    # plot and indicate all peaks
    fig: Figure = plt.figure(figsize=(11.7, 8.3), dpi=100)
    ax: Axes = fig.subplots()
    fig.subplots_adjust(left=0.07, right=0.93, bottom=0.1,
                        top=0.90, wspace=0.05, hspace=0.05)

    np_peaks: npt.NDArray = np.array(
        [], dtype=[('cID', 'i8'), ('Start', 'f8'), ('End', 'f8'), ('Area', 'f8')])
    if args.manual:
        np_peaks = np.genfromtxt(peaks_fileName, names=True,
                                 dtype=['i8', 'f8', 'f8', 'f8'])

    # Automatically Intergate the peaks if use the original data of fid file
    if args.auto:
        peak_list: list = []
        sorted_cID_peaks: npt.NDArray = np.sort(new_peaks, order='cID')

        for cID in set(sorted_cID_peaks['cID']):
            intp_cID: npt.NDArray[np.intp] = (
                np.argwhere(sorted_cID_peaks['cID'] == cID))
            l_Axis = sorted_cID_peaks[intp_cID.min()]['X_AXIS']
            r_Axis = sorted_cID_peaks[intp_cID.max()]['X_AXIS']
            l_LW = sorted_cID_peaks[intp_cID.min()]['X_LW']
            r_LW = sorted_cID_peaks[intp_cID.max()]['X_LW']
            if l_LW <= 1.0:
                l_LW = 1
            if r_LW <= 1.0:
                r_LW = 1
            l_LW_thr = 120/l_LW
            r_LW_thr = 120/r_LW
            l_peak = uc.ppm(l_Axis)+(l_LW/10000)*l_LW_thr
            r_peak = uc.ppm(r_Axis)-(r_LW/10000)*r_LW_thr

            min: int = uc(l_peak, "ppm")
            max: int = uc(r_peak, "ppm")
            if min > max:
                min, max = max, min

            # extract the peak
            peak_list.append((cID, l_peak, r_peak, data[min:max+1].sum()))

        np_peaks = np.array(
            peak_list, dtype=[('cID', 'i8'), ('Start', 'f8'), ('End', 'f8'), ('Area', 'f8')])

    # Intergrate the peaks if manually fixed the plot_1D_peaks.out file
    if args.manual:
        idx: int = 0
        for cID, start, end, _ in np_peaks:  # type: ignore
            min: int = uc(start, "ppm")
            max: int = uc(end, "ppm")
            if min > max:
                min, max = max, min
            # extract the peak
            peak = data[min:max + 1]
            peak_scale = uc.ppm_scale()[min:max + 1]
            # ic(cID, start, end, peak.sum())
            np_peaks['Area'][idx] = peak.sum()  # type: ignore
            idx += 1

    # Delete or Merge peaks if manually read the file
    if args.manual:
        if args.delete:
            for x in args.delete:
                if x in np_peaks['cID']:
                    np_peaks = np_peaks[np_peaks['cID'] != x]
                else:
                    print("  Delete element is wrong cID")
                    print("  Exit and Close the program !!!")
                    exit(1)

        if args.merge:
            min_cID = np.array(args.merge).min()
            start, end = -99999, 99999
            for x in sorted(args.merge):
                if x in np_peaks['cID']:
                    args_x: npt.NDArray[np.intp] = np.argwhere(
                        np_peaks['cID'] == x)
                    # ic(args_x)
                    # ic(np_peaks[args_x[0]]['Start'])
                    if start < np_peaks[args_x[0]]['Start']:
                        start = np_peaks[args_x[0]]['Start'].max()

                    # ic(np_peaks[args_x[0]]['End'])
                    if end > np_peaks[args_x[0]]['End']:
                        end = np_peaks[args_x[0]]['End'].min()
                    # ic(np_peaks[args_x][0])
                else:
                    print("your merge element is wrong cID")
                    print("  Exit and Close the program !!!")
                    exit(1)
            # ic(min_cID, start, end)

            # Remove unnecessary entry
            args_x = np.argwhere(np_peaks['cID'] == min_cID)
            for x in (args.merge):
                if (x in np_peaks['cID']):
                    np_peaks = np_peaks[np_peaks['cID'] != x]

            # assign the new cID of data

            min: int = uc(start, "ppm")  # type: ignore
            max: int = uc(end, "ppm")  # type: ignore
            np_peaks = np.insert(
                np_peaks, args_x[0], (min_cID, start, end, data[min:max+1].sum()))
            # ic(np_peaks)

        if args.cut:

            if args.cut in np_peaks['cID']:
                # use ng.peakpick.pick from y_heighest 0.99 to down two different peaks
                args_x = np.argwhere(np_peaks['cID'] == args.cut)
                l_peaks = np_peaks[args_x][0]['Start']
                r_peaks = np_peaks[args_x][0]['End']
                # ic(l_peaks, r_peaks)
                min: int = uc(l_peaks[0], "ppm")
                max: int = uc(r_peaks[0], "ppm")
                if min > max:
                    min, max = max, min
                # y_cut_heighest = data[min:max+1].max()*0.99
                cut_peaks = ng.peakpick.pick(
                    data=data[min:max+1], pthres=threshold, algorithm="downward")
                sorted_cut_peaks = np.sort(cut_peaks, order='VOL')
                # factor = 1
                # ic(cut_peaks)
                # ic(sorted_cut_peaks)
                end = int(sorted_cut_peaks['X_AXIS'][-1]+min)
                start = int(sorted_cut_peaks['X_AXIS'][-2]+min)
                # ic(data[start:end+1])
                cut_argmin = np.argmin(data[start:end+1])
                cut_center = uc.ppm(start+cut_argmin)
                # ic(cut_center)

                # remove the old entry and add two additional entry
                np_peaks = np_peaks[np_peaks['cID'] != args.cut]
                # ic(start+cut_argmin, min)
                # ic(max, start+cut_argmin)
                np_peaks = np.insert(
                    np_peaks, args_x[0], (args.cut, cut_center, r_peaks, data[min:start+cut_argmin].sum()))
                np_peaks = np.insert(
                    np_peaks, args_x[0], (args.cut, l_peaks, cut_center, data[start+cut_argmin:max].sum()))
            else:
                print("  Merge element is wrong cID")
                print("  Exit and Close the program !!!")
                ic()
                exit(1)

            # Draw the intergral lines
    if args.auto or args.manual:
        for cID, start, end, _ in np_peaks:  # type: ignore
            min: int = uc(start, "ppm")
            max: int = uc(end, "ppm")
            if min > max:
                min, max = max, min

            # extract the peak
            peak = data[min:max + 1]
            peak_scale = uc.ppm_scale()[min:max + 1]
            # plot the integration lines, limits and name of peaks
            ax.plot(peak_scale, peak.cumsum() / 100./4 + peak.max()*0.8, 'g-')
            # ic(cID, peak.sum())
            # ax.plot(peak_scale, [0] * len(peak_scale), 'r-')
            ax.text(peak_scale[0], 0.5 * peak.sum() / 100./4 + peak.max()*0.8, cID,
                    fontsize=8)

    # Saved the plot_1D.peaks
    if args.save:
        np.savetxt(peaks_fileName, np_peaks, fmt=['%10i', '%12.6f', '%12.6f', '%15.1f'], header="     cID        Start          End        Area")  # type: ignore  # nopep8

    # add markers for peak positions. It is not necessary
    if not args.auto and not args.manual:
        for peak, cID, _, _ in new_peaks:
            height = data[int(peak)]
            ppm_peak = uc.ppm(peak)
            if ppm_peak < args.end and ppm_peak > args.start:
                ax.scatter(ppm_peak, height, marker="o", color="r", s=100, alpha=0.5)  # type: ignore # nopep8
                # ax.text(ppm, height + threshold*5, f"{contour_heights[n]:12.3f}" , ha="center", va="center",rotation=90)
                ax.text(ppm_peak, height*1.05, str(cID),
                        ha="center", va="center")
            # print(f"{cID:6d} {ppm:>15.5f}")

    # draw the spectra
    plt.plot(uc.ppm_scale(), data, 'b', linewidth=1)

    # draw the threshold line
    if not args.manual:
        plt.hlines(threshold, args.end, args.start, linestyles="--")  # type: ignore # nopep8
    # plt.hlines(threshold, args.end, args.start,
    #           linestyles="dashdot", linewidth=0.5)

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
    fig.suptitle(args.path, fontsize=12, y=0.98)
    fig.text(0.5, 0.04, "$\\delta$ / ppm", ha="center", fontsize=12)

    from censo_ext.Tools.utility import save_figure
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(pipe_fid_filename)  # type: ignore
    save_figure()

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
