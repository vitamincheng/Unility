#!/usr/bin/env python
import argparse
from icecream import ic
from sys import argv as sysargv
import numpy as np
import numpy.typing as npt
from pathlib import Path

descr = """
________________________________________________________________________________
| dat2BOBYQA.py  
| Usages   : dat2BOBYQA.py <dat file> [options]
| Input    : -i input file [default anmr.dat]
| [options]
| Dir      : -D Directory [default .]
| Start    : -s --start start point in spectra [default -5]
| End      : -e --end end point in spectra [default 15]
| dpi      : --dpi 100 for H 5 for C [defalt 100]
| If your spectra is below 50 ppm and args.end is below 50 ppm, it will 
| automatically set the parameter to start=-20ppm end=250ppm and dpi=500
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-D",
        "--dir",
        dest="dir",
        action="store",
        required=False,
        default=".",
        help="Provide the directory name [default .]",
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="anmr.dat",
        help="Provide input_file name [default anmr.dat]",
    )

    parser.add_argument(

        "-s",
        "--start",
        dest="start",
        action="store",
        required=False,
        type=float,
        default=-5.0,
        help="Provide start point ppm [default -5]",
    )

    parser.add_argument(
        "-e",
        "--end",
        dest="end",
        action="store",
        required=False,
        type=float,
        default=15.0,
        help="Provide end point ppm [default 15]",
    )

    parser.add_argument(
        "--dpi",
        dest="dpi",
        action="store",
        required=False,
        type=int,
        default=100,
        help="dpi 100(for H), 5(for C) [default 100]",
    )

    parser.add_argument(
        "-t",
        "--thr",
        dest="thr",
        action="store",
        required=False,
        type=int,
        default=10,
        help="threshold of baseline [default 10]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print(f"    provided arguments: {" ".join(sysargv)}")

    if args.file:
        from censo_ext.Tools.utility import IsExist
        IsExist(args.file)
        from censo_ext.Tools.datfile import CensoDat
        inDat: CensoDat = CensoDat(args.file)
        inDat.method_normalize_dat(
            start=args.start, end=args.end, dpi=args.dpi)
        # inDat.set_fileName(args.out)
        # inDat.method_save_dat()
        from scipy.signal import find_peaks
        x: npt.NDArray[np.float64] = inDat.get_Dat()
        from censo_ext.Tools.spectra import numpy_thr
        thr: float = numpy_thr(x.T[1], args.thr)
        peaks, _ = find_peaks(x.T[1], thr)
        peaks_np: npt.NDArray[np.float64] = x.T[0][peaks]*(-1)
        # ic(peaks_np)

        DirFileName: Path = args.dir / Path("Average/NMR/orcaS.out")
        from censo_ext.Tools.utility import IsExists_DirFileName
        IsExists_DirFileName(DirFileName)
        np_data = np.genfromtxt(DirFileName)
        # ic(np_data)
        from censo_ext.Tools.spectra import find_nearest

        app_list: list[float] = []
        peaks_list = [x[1] for x in np_data]
        for x in peaks_list:
            value, _ = find_nearest(list(peaks_np), x)
            app_list.append(value)
        # ic(app_list)


if __name__ == "__main__":
    main()

#
#
#
