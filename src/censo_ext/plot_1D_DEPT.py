#! /usr/bin/env python
from pathlib import Path
from icecream import ic
import nmrglue as ng
from nmrglue.fileio.fileiobase import unit_conversion
import numpy as np
import numpy.typing as npt
import argparse
import sys
from censo_ext.Tools.utility import delete_all_files, print_arguments


descr = """
________________________________________________________________________________
| For Plot 1D sepctra in experiments using nmrglue module
| Usage: plot_1D_DEPT.py <geometry> [options]
| [Options]
| Dir      : -d Directory of path of 13C /DEPT_135 and DEPT_90 [default None]
|          : -start start point of chemical shift [default from data]
|          : -end   end point of chemical shift [default from data]
| Save     : --save saved the report of carbon [default false]
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

    return parser.parse_args()


# global variable
peaks_fileName = "plot_1D_DEPT.peaks"
fid_fileName = ".1d_pipe.fid"


def Channel(args: argparse.Namespace, path: dict[Path, Path], channel: Path, thr: float, phase: float = 1.0) -> dict[int, float]:
    dic, data = ng.bruker.read_pdata(str(path[channel]))
    udic = ng.bruker.guess_udic(dic, data)

    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    ng.pipe.write(fid_fileName, *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(fid_fileName)
    data = data.real*phase  # type: ignore
    uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)

    # end ---------+--------- start
    # args.end                args.start
    # ppm_1h_0                ppm_1h_1

    ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    ppm = np.linspace(ppm_1h_0, ppm_1h_1, data.shape[0])
    if not args.start or not args.end:
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
    Result: dict = dict()

    for n, peak in enumerate(peaks):
        ppm: float = uc_1h.ppm(peak)
        Result[n+1] = ppm
    return Result


def Compare_two_dict(CH1: dict[int, float], CH2: dict[int, float], StAtoms: dict[int, int], Label: int) -> None:
    from censo_ext.Tools.spectra import find_nearest
    for x in CH2.values():
        nearest_peak, idx0 = find_nearest(list(CH1.values()), x)
        if (nearest_peak-x) < 0.05:
            StAtoms[idx0+1] = Label
        else:
            print("some peaks is more than 0.02 ppm")
            print("  Exit and Close the program !!!")
            ic()
            exit(0)


def print_outcome(StAtoms, uc_1h, peaks) -> None:
    print("#   ID             ppm    nHydrogens")
    for n, peak in enumerate(peaks):
        ppm = uc_1h.ppm(peak)
        print(f"{n+1:6d} {ppm:>15.5f}        {StAtoms[n+1]:>3d}")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    # read in the Bruker data
    if args == argparse.Namespace():
        args = cml()
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

    path[ch1] = directory / ch1 / Path("pdata/1")
    path[ch2] = directory / ch2 / Path("pdata/1")
    path[ch3] = directory / ch3 / Path("pdata/1")

    channel: Path = ch1
    Result_ch1: dict[int, float] = Channel(args, path, channel, thr[channel])
    StAtoms: dict[int, int] = {key: -1 for key, value in Result_ch1.items()}

    channel: Path = ch2
    Result_ch2: dict[int, float] = Channel(args, path, channel, thr[channel])
    print("DEPT90             ppm")
    for idx1, ch2_ppm in Result_ch2.items():
        print(f"{idx1:6d} {ch2_ppm:>15.5f}")

    channel: Path = ch3
    Result_ch3: dict[int, float] = Channel(args, path, channel, thr[channel])
    Compare_two_dict(Result_ch1, Result_ch3, StAtoms, Label=3)
    Compare_two_dict(Result_ch1, Result_ch2, StAtoms, Label=1)
    print("DEPT135(up)        ppm")
    for idx1, ch3_ppm in Result_ch3.items():
        print(f"{idx1:6d} {ch3_ppm:>15.5f}")

    channel: Path = ch3
    Result_ch3_180: dict[int, float] = Channel(args, path, channel, thr_ch3_180, phase=-1.0)  # nopep8
    Compare_two_dict(Result_ch1, Result_ch3_180, StAtoms, Label=2)
    print("DEPT135(down)      ppm")
    for idx1, ch3_180_ppm in Result_ch3_180.items():
        print(f"{idx1:6d} {ch3_180_ppm:>15.5f}")
    for key, value in StAtoms.items():
        if value == -1:
            StAtoms[key] = 0

    channel = ch1
    dic: dict
    data: npt.NDArray
    dic, data = ng.bruker.read_pdata(str(path[channel]))
    udic: dict = ng.bruker.guess_udic(dic, data)

    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    ng.pipe.write(fid_fileName, *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(fid_fileName)
    data = data.real  # type: ignore
    uc_1h: unit_conversion = ng.pipe.make_uc(dic, data)

    # end ---------+--------- start
    # args.end                args.start
    # ppm_1h_0                ppm_1h_1

    # ppm_1h_0, ppm_1h_1 = uc_1h.ppm_limits()
    # ppm: npt.NDArray[np.float64] = np.linspace(
    #    ppm_1h_0, ppm_1h_1, data.shape[0])

    if not args.start or not args.end:
        args.end, args.start = uc_1h.ppm_limits()

    from censo_ext.Tools.spectra import numpy_thr_mean_3
    threshold: float = numpy_thr_mean_3(
        data.astype(np.float64))*thr[channel]

    # detect all peaks with a threshold
    from scipy.signal import find_peaks
    y_heighest = max(data)
    threshold += float(y_heighest) * 0.01
    peaks, _ = find_peaks(data, height=threshold, width=1)

    # print the final data
    print_outcome(StAtoms, uc_1h, peaks)

    # save to file
    if args.save:
        with open(peaks_fileName, "w") as f:
            sys.stdout = f
            print_outcome(StAtoms, uc_1h, peaks)
        sys.stdout = sys.__stdout__

    delete_all_files(fid_fileName)


if __name__ == "__main__":
    main()
