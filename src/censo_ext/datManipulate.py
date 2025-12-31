#!/usr/bin/env python
from pathlib import Path
from censo_ext.Tools.utility import print_arguments
import argparse
descr = """
________________________________________________________________________________
| For Generation of xyz molecule 
| Usages   : xyz.py <geometry> [options]
| [options]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="descr",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="anmr.dat",
        help="Provide one input xyz file [default anmr.dat]",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.dat",
        help="Provide output dat file name [default output.dat]",
    )
    parser.add_argument(
        "--shift",
        dest="shift",
        action="store",
        required=False,
        type=float,
        default=0.0,
        help="Provide the shift ppm [default 0.0]",
    )
    parser.add_argument(
        "--intensit",
        dest="intensit",
        action="store",
        required=False,
        type=float,
        default=1.0,
        help="Provide the shift ppm [default 1.0]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile: Path = Path(args.file)
    outFile: Path = Path(args.out)
    if inFile:
        from censo_ext.Tools.utility import IsExist
        IsExist(inFile)
        from censo_ext.Tools.datfile import CensoDat
        inDat: CensoDat = CensoDat(inFile)
        inDat.ppm_shift_Dat(ppm_shift=args.shift)
        inDat.Intensit_Dat(Intensit=args.intensit)
        inDat.set_fileName(outFile)
        inDat.method_save_dat()


if __name__ == "__main__":
    main()
