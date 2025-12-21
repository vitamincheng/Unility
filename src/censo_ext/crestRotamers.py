#!/usr/bin/env python
from pathlib import Path
from censo_ext.Tools.utility import IsExist_bool, print_arguments
import argparse
import numpy as np
import numpy.typing as npt

from censo_ext.Tools.xyzfile import GeometryXYZs
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
        default="crest_rotamers.xyz",
        help="Provide one input xyz file [default crest_rotamers.xyz]",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile: Path = Path(args.file)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    fileName = Path("cre_members")

    if IsExist_bool(fileName):
        np_inData: npt.NDArray[np.int64] = np.genfromtxt(
            fileName, skip_header=1, dtype=int)

        for _, start, end in np_inData:
            import copy
            outFile: GeometryXYZs = copy.deepcopy(xyzFile)
            outFile.method_xyzExtract([*range(start-1, end)])
            outFile.set_filename(f"{start}_{end}.xyz")
            outFile.method_save_xyz([])
            print(f"Data saved to : {start}_{end}.xyz")
    else:
        print(f"  Your file {fileName} is not Exist !!! ")
        print("  Close and Exit the program !!!")
        exit(0)


if __name__ == "__main__":
    main()
