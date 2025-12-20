#!/usr/bin/env python
import argparse
# import os
from scipy.spatial.transform import Rotation as R
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.factor import idx_3atom_opt
from censo_ext.Tools.utility import AtomID, print_arguments
from pathlib import Path
from censo_ext.Tools.xyzfile import GeometryXYZs

descr = """
________________________________________________________________________________
| For Return to origin and lay on the xz plane
| Usage    : xyzReturnOandZ.py <geometry> [options]                  
| Input    : -i xyz file [default traj.xyz]
| Output   : -o Saved xyz file [default output.xyz] 
| [Options]
| Atom     : -a or --atom [1 2 3] idx of atom 
|              1 : Fixed atom and return origin
|              2 : Rotation atom and z axis
|              3 : Rotation atom and lay on xz plane
| Automatic: --auto Automatically search for the minimum deviation setting 
|                   origin of all atoms [default False] 
| Replace  : -r Replace the input file [default False]
| Print    : -p Print the final data on screen [default False]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description=f"{descr}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

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
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Provide one output xyz file [default output.xyz]",
    )

    parser.add_argument(
        "-a",
        "--atom",
        dest="atom",
        action="store",
        type=int,
        nargs=3,
        help="Provide three idx of atom's nubmers [1(origin) 2(z axis) 3(xz plane)]",
    )

    parser.add_argument(
        "-r",
        "--replace",
        dest="replace",
        action="store_true",
        help="Replace the original input file [default False]",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        help="Print the final data on stdout (on screen) [default False]",
    )

    parser.add_argument(
        "--auto",
        dest="auto",
        action="store_true",
        help="Automatically search for the minimum deviation setting origin of all atoms [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    from censo_ext.Tools.utility import IsExist
    inFile = Path(args.file)
    if args.replace:
        outFile: Path = inFile
    else:
        outFile = Path(args.out)

    IsExist(args.file)

    if not args.atom and not args.auto:
        raise ValueError(" No any sepific atom in your provided arguments ")

    p_idx1: AtomID
    q_idx1: AtomID
    r_idx1: AtomID

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    if not args.atom and args.auto:
        if len(xyzFile.Sts) == 1 or 0:
            print(
                f"  {xyzFile.get_fileName()}, the numbers of the structures is {len(xyzFile)}")
            print("  The numbers of the structures of your xyzFile is one or zero")
            print("  Not Need to find the 3 numbers of opt atomic index !!!")
            exit(0)
        print("\n Automated to set the 3 atoms to return origin and lay on XZ plane")
        print(" First FactorAnalysis.py will executive and second continue the RetrunOandZ.py ")
        p_idx1, q_idx1, r_idx1 = idx_3atom_opt(xyzFile)
    else:
        p_idx1, q_idx1, r_idx1 = args.atom

    # Process xyz file
    for St in xyzFile.Sts:

        dxyz: npt.NDArray[np.float64] = St.coord[p_idx1-1].copy()
        St.coord -= dxyz  # type: ignore
        z_axis = (0, 0, np.sqrt(
            np.sum(np.square(St.coord[q_idx1-1]))))

        rotation_axis = St.coord[q_idx1-1] + z_axis

        Normalized_RotationAxis: npt.NDArray[np.float64] = np.array([
            0, 1, 0]) if np.linalg.norm(rotation_axis) == 0 else rotation_axis / np.linalg.norm(rotation_axis)

        R_pq = R.from_rotvec(np.pi*Normalized_RotationAxis)
        St.coord = R_pq.apply(St.coord)  # type: ignore

        Angle_qr = np.angle(
            complex(St.coord[r_idx1-1][0], complex(St.coord[r_idx1-1][1])))
        R_qr = R.from_euler('z', -Angle_qr)
        St.coord = R_qr.apply(St.coord)  # type: ignore

    # Save or print result
    if args.print:
        xyzFile.method_print([])
    else:
        fileName: Path = inFile if args.replace else outFile
        print(f"    Saved to {fileName}")
        xyzFile.set_filename(fileName)
        xyzFile.method_save_xyz([])


if __name__ == "__main__":
    main()

    # test
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz -a 30 45 47
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz --auto
