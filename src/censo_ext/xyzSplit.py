#!/usr/bin/env python3
from scipy.spatial.transform import Rotation as R
import argparse
import os
import numpy as np
from censo_ext.Tools.xyzfile import GeometryXYZs
from sys import argv as sysargv
descr = """
________________________________________________________________________________
|                                          [07.19.2023] vitamin.cheng@gmail.com
| For search the confomrers from various angles of cleavage specifying two atoms                        
| Usage    : xyzSplit.py [options]                  
| [Options]
| Input    : -i Read xyz file [default traj.xyz]
| Atom     : -a or --atom [1 2] idx of atom's number  
|              1 : Fixed atom
|              2 : Rotation axis atom (360 degrees) 
| nCut     : -c or cut Number of cut to make 360 degrees around the roation axis 
| Output   : -o Save xyz file [default output.xyz] 
| Print    : -p Print output to screen
| Packages : Tools 
| Module   : xyzfile.py / topo.py / unility.py
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
        help="Provide one output xyz file [defalut output.xyz]",
    )

    parser.add_argument(
        "-a",
        "--atom",
        dest="atoms",
        action="store",
        type=int,
        nargs=2,
        default=None,
        required=False,
        metavar=('FIXED', 'ROTATION'),
        help="two atom indics: first is fixed, second is rotation axis"
    )

    parser.add_argument(
        "-c",
        "--cut",
        dest="cuts",
        action="store",
        type=int,
        default=None,
        required=False,
        help="Number of cuts to make in 360 degrees around the rotation axis",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        default=False,
        help="Print output to screen",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)

    from censo_ext.Tools.utility import is_exist
    from pathlib import Path
    if is_exist(Path(args.file)):
        pass

    if not args.print:
        print(descr)  # Program description
        print("    provided arguments: {}".format(" ".join(sysargv)))

    if args.cuts == None or args.atoms == None:
        print(" Please input your atoms that you want to split ")
        print(" Exit to this program !!!")
        exit(0)

    from censo_ext.Tools.utility import delete_file_bool
    Delete_work: bool = False
    if not args.print:
        Delete_work = delete_file_bool(args.out)

    idx1_p, idx1_q = args.atoms[0], args.atoms[1]
    nCutters: int = args.cuts

    x: dict = {"file": args.file, "bond_broken": [
        idx1_q, idx1_p], "print": False, "debug": False}
    from censo_ext.Tools.topo import Topo
    Sts_topo: Topo = Topo(x["file"])
    idx0_list: list[int] = [
        x-1 for x in Sts_topo.method_broken_bond_H(argparse.Namespace(**x))]

    infile: GeometryXYZs = GeometryXYZs(args.file)
    infile.method_read_xyz()

    for idx0_st, St in enumerate(infile.Sts):

        dxyz: np.ndarray = St.coord[idx1_p-1]
        inital: list[np.ndarray] = St.coord[:].copy()

        for nCutter in range(nCutters):

            St.coord[:] = inital.copy()
            St.coord[:] -= dxyz

            rotation_axis = St.coord[idx1_q-1]
            rotation_vector = rotation_axis/np.linalg.norm(rotation_axis)

            r_pq = R.from_rotvec(2*np.pi*(nCutter/nCutters)*rotation_vector)

            for idx in idx0_list:
                St.coord[idx] = r_pq.apply(St.coord[idx])

            St.coord[:] += dxyz

            if args.print:
                infile.method_print([idx0_st+1])
            else:
                # ic(idx0_st,nCutter)
                infile.set_filename(args.out)
                infile.method_save_xyz_append([idx0_st+1])

    if not args.print:
        if Delete_work == True:
            print(f"   {args.out} is here, it will be removed.")
            print(f"    Overwrite the file : {args.out}")
        else:
            print(f"    Create a new file  : {args.out}")


if __name__ == "__main__":
    main()

    # For example
    #   python3 xyzSplit.py -i ../tests/crest_conformers1.xyz -a 52 55 -c 3
    #   python3 xyzSplit.py -i ../tests/crest_conformers1.xyz -a 52 55 -c 3 -p
