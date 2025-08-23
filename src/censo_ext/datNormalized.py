#!/usr/bin/env python
# from icecream import ic
import argparse
from sys import argv as sysargv
# from scipy.signal import argrelmax, argrelmin

descr = """
________________________________________________________________________________
|                                          [09.16.2024] vitamin.cheng@gmail.com
| dat_normalized.py  
| Usages  : dat_normalized.py [options]
| [options]
| Input    : -i input file [default anmrh.dat]
| Output   : -o output file [default output.dat]
| Start    : -s --start start point in spectra [default -5]
| End      : -e --end end point in spectra [default 15]
| dpi      : --dpi 10000 for H 500 for C [defalt 10000]
| Package  : Tools 
| Module   : anmrfile.py / unility.py
| If your spectra is below 50 ppm and args.end is below 50 ppm, it will 
| automatically set the parameter to start=-20ppm end=250ppm and dpi=500
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="anmrh.dat",
        help="Provide input_file name [default anmrh.dat]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.dat",
        help="Provide output_file name [default output.dat]",
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
        default=10000,
        help="dpi 10000(for H), 500(for C) [default 10000]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.file:
        from censo_ext.Tools.utility import IsExist
        IsExist(args.file)
        from censo_ext.Tools.datfile import CensoDat
        inDat: CensoDat = CensoDat(args.file)
        inDat.method_normalize_dat(
            start=args.start, end=args.end, dpi=args.dpi)
        inDat.set_fileName(args.out)
        inDat.method_save_dat()


if __name__ == "__main__":
    main()
