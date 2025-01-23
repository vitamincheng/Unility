#!/usr/bin/env python3
from scipy.spatial.transform import Rotation as R
import argparse
import numpy as np
from Tools.xyzfile import ClassGeometryXYZs
from icecream import ic
from sys import argv as sysargv
descr = """
________________________________________________________________________________
|                                          [07.19.2023] vitamin.cheng@gmail.com
| For search the confomrers from various angles of cleavage specifying two atoms                        
| Usage    : xyzSplit.py [options]                  
| [Options]
| Input    : -i xyz file [default traj.xyz]
| Atom     : -a or --atom [1 2] idx of atom's number  
|              1 : Fixed atom
|              2 : Rotation atom (360 degrees) 
| nCut     : -c or cut numbers cutting numbers of atom in 360 degrees 
| Output   : -o Saved to xyz file [default output.xyz] 
| Print    : -p Print the final data on screen
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
        dest="atom",
        action="store",
        type=int,
        nargs=2,
        default=None,
        required=False,
        help="[1,2] 1: fixed atom 2: rotation atom. Provide two idx of atom's nubmers"
    )

    parser.add_argument(
        "-c",
        "--cut",
        dest="cut",
        action="store",
        type=int,
        default=None,
        required=False,
        help="Provide the total cutting numbers in your atom's in 360 degrees",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        default=False,
        help="Print the final data on screen (stdout)",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    if not args.print:
        print(descr)  # Program description
        print("    provided arguments: {}".format(" ".join(sysargv)))
    # ic(args)
    if args.cut == None or args.atom == None:
        print(" Please input your atoms that you want to split ")
        print(" Exit to this program !!!")
        import os
        ic()
        os._exit(0)

    from Tools.utility import DeleteFileBool
    Delete_work: bool = False
    if not args.print:
        Delete_work = DeleteFileBool(args.out)

    idx1_p, idx1_q = args.atom[0], args.atom[1]
    nCutters: int = args.cut

    x: dict = {"file": args.file, "bond_broken": [
        idx1_q, idx1_p], "print": False, "debug": False}
    from Tools.topo import topo
    Sts_topo: topo = topo(x["file"])
    idx0_list: list[int] = [
        x-1 for x in Sts_topo.Broken_bond_H(argparse.Namespace(**x))]
    # ic(idx1_list)
    infile: ClassGeometryXYZs = ClassGeometryXYZs(args.file)
    infile.method_read_xyz()
    # outfile = zyzfile.ClassGeometryXYZs(args.out)
    # ic(infile.get_nSt(),nCutters)
    for idx0_st in range(len(infile)):

        dxyz: np.ndarray = infile.Sts[idx0_st].coord[idx1_p-1].copy()
        inital: list[np.ndarray] = infile.Sts[idx0_st].coord[:].copy()

        for nCutter in range(nCutters):

            infile.Sts[idx0_st].coord[:] = inital.copy()
            infile.Sts[idx0_st].coord[:] -= dxyz

            pq_vector = infile.Sts[idx0_st].coord[idx1_q-1]
            n_pq_vector = pq_vector/np.linalg.norm(pq_vector)

            r_pq = R.from_rotvec(2*np.pi*(nCutter/nCutters)*n_pq_vector)

            for idx in idx0_list:
                infile.Sts[idx0_st].coord[idx] = r_pq.apply(
                    infile.Sts[idx0_st].coord[idx])

            infile.Sts[idx0_st].coord[:] += dxyz

            if args.print:
                infile.method_print([idx0_st+1])
            else:
                # ic(idx0_st,nCutter)
                infile.set_filename(args.out)
                infile.method_save_xyz_append([idx0_st+1])

    if not args.print:
        if Delete_work == True:
            print("   ", args.out, " is here, it will be removed.")
            print("    Overwrite the file : ", args.out)
        else:
            print("    Create a new file : ",  args.out)


if __name__ == "__main__":
    main()

    # For example
    #   python3 xyzSplit.py -i ../tests/crest_conformers1.xyz -a 52 55 -c 3
    #   python3 xyzSplit.py -i ../tests/crest_conformers1.xyz -a 52 55 -c 3 -p
