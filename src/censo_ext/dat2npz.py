#!/usr/bin/env python
import argparse
import numpy as np
from censo_ext.Tools.utility import IsExists_DirFileName
descr = """
________________________________________________________________________________
| For Transform from dat to npz file and reverse 
| Usages   : dat2npz.py <geometry> [options]
| [options]
| input    : -i input dat or npz file 
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=True,
        help="Provide one input dat or npz file",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    if args.file:
        path, file = IsExists_DirFileName(args.file)
        file_split = file.split(".")
        file_ext = file_split[1]
        fileName = file_split[0]

        if file_ext == "dat":
            in_Data = np.genfromtxt(args.file)
            np.savez_compressed(fileName+".npz", in_Data)
            print(f" the spectra is saved to : {fileName+'.npz'}")
        elif file_ext == "npz":
            in_Data = np.load(args.file)
            np.savetxt(fileName+".dat", in_Data["arr_0"], fmt='%2.5f %12.5e')
            print(f" the spectra is saved to : {fileName+'.dat'}")
        else:
            print("  The file is not supported for this system !!!")
            print("  Exit and Close the program !!!")
            exit(0)


if __name__ == "__main__":
    main()
