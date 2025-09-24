#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.utility import IsExists_DirFileName
from censo_ext.Tools.utility import delete_all_files
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
        file_split: list[str] = file.split(".")
        file_ext: str = file_split[1]
        fileName: str = file_split[0]

        if file_ext == "dat":
            in_Data: npt.NDArray[np.float64] = np.genfromtxt(args.file)
            np.savez_compressed(fileName + ".npz", in_Data)
            print(f" the spectra is saved to : {fileName + '.npz'}")
            delete_all_files(args.file)
        elif file_ext == "npz":
            in_Data: npt.NDArray[np.float64] = np.load(args.file)["arr_0"]
            np.savetxt(fileName + ".dat", in_Data, fmt='%2.5f %12.5e')
            print(f" the spectra is saved to : {fileName + '.dat'}")
            delete_all_files(args.file)
        else:
            print("  The file is not supported for this system !!!")
            print("  Exit and Close the program !!!")
            exit(0)


if __name__ == "__main__":
    main()
