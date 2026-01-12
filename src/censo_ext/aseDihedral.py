#!/usr/bin/env python
from censo_ext.Tools.utility import AtomID
from pathlib import Path
from censo_ext.Tools.utility import print_arguments
import argparse
import matplotlib.pyplot as plt
from censo_ext.Tools.ml4nmr import ase_get_dihedral
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
        default="traj.xyz",
        help="Provide one input xyz file [default traj.xyz]",
    )
    parser.add_argument(
        "-a",
        "--atom",
        dest="atom",
        action="store",
        type=int,
        nargs=4,
        help="Provide three idx1 of atom's nubmers of Dihedral",
    )
    args: argparse.Namespace = parser.parse_args()
    return args


type cell_4AtomIDs = tuple[AtomID, AtomID, AtomID, AtomID]


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile = Path(args.file)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    degree: list[float] = []
    atoms: list[int] = args.atom-1

    for idx1 in range(1, len(xyzFile)+1):
        in_cell: cell_4AtomIDs = (AtomID(atoms[0]), AtomID(
            atoms[1]), AtomID(atoms[2]), AtomID(atoms[3]))
        res: float = ase_get_dihedral(
            xyzFile=xyzFile, idx1=idx1, in_cell=in_cell)
        degree.append(float(res))
    # print(degree)
    import numpy as np
    print(len(degree))
    print(np.average(np.array(degree)))

    plt.hist(degree, bins=24, density=False, color='blue',
             edgecolor='black')  # density=False plots counts

    # 3. Add labels and a title
    plt.ylabel('Frequency')
    plt.xlabel('Data Values')
    plt.title('Distribution of Data (Histogram)')

    # 4. Display the plot
    plt.show()


if __name__ == "__main__":
    main()
