#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
from pathlib import Path
from censo_ext.Tools.utility import print_arguments

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


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        # description=f"{descr}",
        # formatter_class=argparse.RawDescriptionHelpFormatter,
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

    return parser.parse_args()


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    if args.file:
        from censo_ext.Tools.utility import IsExist
        IsExist(args.file)

        from censo_ext.Tools.datfile import CensoDat
        inDat: CensoDat = CensoDat(args.file)
        inDat.method_normalize_dat(
            start=args.start, end=args.end, dpi=args.dpi)

        from scipy.signal import find_peaks
        x: npt.NDArray[np.float64] = inDat.get_Dat()

        from censo_ext.Tools.spectra import numpy_thr
        thr: float = numpy_thr(x.T[1], args.thr)
        peaks, _ = find_peaks(x.T[1], thr)
        peaks_np: npt.NDArray[np.float64] = x.T[0][peaks]*(-1)

        DirFileName: Path = args.dir / Path("Average/NMR/orcaS.out")
        from censo_ext.Tools.utility import IsExists_DirFileName
        IsExists_DirFileName(DirFileName)
        np_data: npt.NDArray = np.genfromtxt(DirFileName)

        from censo_ext.Tools.spectra import find_nearest

        app_list: list[float] = []
        peaks: list = [x[1] for x in np_data]
        for x in peaks:
            value, _ = find_nearest(peaks_np, x)
            app_list.append(value)


if __name__ == "__main__":
    main()

#
#
#
