#!/usr/bin/env python3
import argparse
import os
from sys import argv as sysargv
from typing import Any
from scipy.spatial.transform import Rotation as R
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.xyzfile import GeometryXYZs
from icecream import ic
from pathlib import Path

descr = """
________________________________________________________________________________
|                                          [10.05.2024] vitamin.cheng@gmail.com
| For Return to origin and lay on the xz plane
| Usage: xyzReturnOandZ.py <geometry> [options]                  
| [Options]
| Input    : -i xyz file [default traj.xyz]
| Atom     : -a or --atom [1 2 3] idx of atom  
|              1 : Fixed atom and return origin
|              2 : Rotation atom and z axis
|              3 : Rotation atom and lay on xz plane
| Automatic: --auto Automatically search for the minimum deviation setting 
|                   origin of all atoms [default False] 
| Output   : -o Saved xyz file [default output.xyz] 
| Replace  : -r Replace the input file
| Print    : -p Print the final data on screen 
| Packages : Tools  
| Module   : xyzfile.py / topo.py / factor.py 
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        #        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )  # argparse.RawDescriptionHelpFormatter) #,

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
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
        required=False,
        help="Provide three idx of atom's nubmers [1(origin) 2(z axis) 3(xz plane)]",
    )

    parser.add_argument(
        "-r",
        "--replace",
        dest="replace",
        action="store_true",
        default=False,
        help="Replace the original input file",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        default=False,
        help="Print the final data on stdout (on screen)",
    )

    parser.add_argument(
        "--auto",
        dest="auto",
        action="store_true",
        default=False,
        help="Automatically search for the minimum deviation setting origin of all atoms [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def idx_3atom_opt(args) -> tuple[int, int, int]:
    from censo_ext.Tools.factor import method_factor_analysis
    args_x: dict = {"file": args.file,
                    "factor": 0.5, "debug": False, "opt": False}
    minor_list: list[int]
    TableSTD: dict[int, float]
    minor_list, TableSTD = method_factor_analysis(
        args=argparse.Namespace(**args_x))
    idx1_Low_Factor: list[int] = list(minor_list)
    idx1_Atom_list: list[int] = list(TableSTD.keys())
    AtomSTD: list[float] = list(map(float, TableSTD.values()))

    idx1_Bonding: list[list[int]] = []
    for idx, x in enumerate(idx1_Low_Factor):
        from censo_ext.Tools.topo import Topo
        args_x: dict = {"file": args.file, "bonding": x,
                        "print": False, "debug": False}
        Sts_topo: Topo = Topo(args_x["file"])
        idx1_Bonding.append(Sts_topo.method_bonding(
            args=argparse.Namespace(**args_x)))

    idx1_3atom: list[list[int]] = []
    for idx, x in enumerate(idx1_Low_Factor):
        # total numbers >=3 or >2 (one of total numbers is )
        if len(idx1_Bonding[idx]) > 1:
            tmp: list[int] = []
            tmp.append(x)
            for y in idx1_Bonding[idx]:
                tmp.append(y)
            idx1_3atom.append(tmp)

    from itertools import combinations
    idx1_Combine_3atom: list[tuple[int, int, int]] = []
    for idx, x in enumerate(idx1_3atom):
        return_list = list(combinations(x, 3))
        for y in (return_list):
            idx1_Combine_3atom.append(y)

    idx_min_Total_Dev = 0
    min_Total_Dev = 100
    for idx, x in enumerate(idx1_Combine_3atom):
        Total_Dev_Atom = 0
        for y in x:
            Total_Dev_Atom += (AtomSTD[idx1_Atom_list.index(y)])
        if min_Total_Dev > Total_Dev_Atom:
            min_Total_Dev = Total_Dev_Atom
            idx_min_Total_Dev: int = idx

    print("")
    print(" 3 atom idx of lowest total factor : ",
          idx1_Combine_3atom[idx_min_Total_Dev])
    print("")
    return (idx1_Combine_3atom[idx_min_Total_Dev])


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)

    if not args.print:
        print(descr)  # Program description
        print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.atom == None and args.auto == False:
        print(" No any sepific atom in your provided arguments ")
        ic()
        exit(0)

    p_idx: int
    q_idx: int
    r_idx: int
    if args.atom == None and args.auto:
        print("\n Automated to set the 3 atoms to return origin and lay on XZ plane")
        print(" First FactorAnalysis.py will executive and second continue the RetrunOandZ.py ")
        p_idx, q_idx, r_idx = idx_3atom_opt(args)
    else:
        p_idx, q_idx, r_idx = args.atom

    infile: GeometryXYZs = GeometryXYZs(args.file)
    infile.method_read_xyz()

    # Process xyz file
    for idx_st in range(len(infile)):

        dxyz: npt.NDArray[np.float64] = infile.Sts[idx_st].coord[p_idx-1].copy()
        infile.Sts[idx_st].coord -= dxyz
        import math
        z_axis = (0, 0, math.sqrt(
            np.sum(np.square(infile.Sts[idx_st].coord[q_idx-1]))))

        rotation_axis = infile.Sts[idx_st].coord[q_idx-1]+z_axis
        if np.linalg.norm(rotation_axis) == 0:
            Normalized_rotation_axis: npt.NDArray[np.float64] = np.array([
                                                                         0, 1, 0])
        else:
            Normalized_rotation_axis = rotation_axis / \
                np.linalg.norm(rotation_axis)

        R_pq = R.from_rotvec(np.pi*Normalized_rotation_axis)
        infile.Sts[idx_st].coord = R_pq.apply(
            infile.Sts[idx_st].coord)

        Angle_qr = np.angle(complex(infile.Sts[idx_st].coord[r_idx-1][0], complex(
            infile.Sts[idx_st].coord[r_idx-1][1])))
        R_qr = R.from_euler('z', -Angle_qr)
        infile.Sts[idx_st].coord = R_qr.apply(
            infile.Sts[idx_st].coord)

    # Save or print result
    if args.print:
        infile.method_print([])
    else:
        filename: Path = Path(args.file) if args.replace else Path(args.out)
        print(f"    Saved to {filename}")
        infile.set_filename(filename)
        infile.method_save_xyz([])


if __name__ == "__main__":
    main()

    # test
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz -a 30 45 47
    # python3 xyzReturnOandZ.py -i ../tests/crest_conformers.xyz --auto
