#!/usr/bin/env python
import argparse
# from sys import argv as sysargv
import subprocess
from icecream import ic
from censo_ext.Tools.xyzfile import GeometryXYZs

descr = """
________________________________________________________________________________
|                                          [09.06.2024] vitamin.cheng@gmail.com
| Generation inp file for PES and Create PES of xtb
| Usage : create_PES.py <geometry> [options]
| [Options]
| Input    : -i one xyz file [default: isomers.xyz]
| Atoms    : -c 3 atoms' idx 1(First) 2(Center) 3(Second)
| Distance : -b1 Bond Distance of 1,2 Atoms From 1 to 2 [unit: Å] 
| Distance : -b2 Bond Distance of 3,2 Atoms From 1 to 2 [unit: Å]
| Cut      : -cut cut part of b1 and b2 each [default: 10,10]
| Max      : -max opt maxmium cycles [default: 10]
|______________________________________________________________________________
"""


def cml():
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
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
        help="Provide one xyz file to reorganize the serial numbers [default: traj.xyz]",
    )

    parser.add_argument(
        "-c",
        "--center-atom",
        dest="center_atom",
        action="store",
        nargs=3,
        required=False,
        type=int,
        help="Location of 3 Atoms : 1(First) 2(Center) 3(Second) ",
    )

    parser.add_argument(
        "-b1",
        "--bond1-distance",
        dest="bond1_distance",
        action="store",
        nargs=2,
        required=False,
        type=float,
        help="Distance of 1's Atom with 2's Atom : 1(From) 2(To)",
    )

    parser.add_argument(
        "-b2",
        "--bond2-distance",
        dest="bond2_distance",
        action="store",
        nargs=2,
        required=False,
        type=float,
        help="Distance of 3's Atom with 2's Atom : 1(From) 2(To)",
    )

    parser.add_argument(
        "-cut",
        "--cut-distance",
        dest="cut_distance",
        action="store",
        nargs=2,
        required=False,
        type=int,
        default=[10, 10],
        help="cut part of bond1 and bond 2 distances [default: 10,10]",
    )

    parser.add_argument(
        "-max",
        "--max-cycle",
        dest="max_cycle",
        action="store",
        required=False,
        type=int,
        default=10,
        help="opt maxmium cycles [default: 10]",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = cml()
    ic(args.center_atom)
    ic(args.bond1_distance)
    ic(args.bond2_distance)
    ic(args.cut_distance)
    ic(args.max_cycle)
    force_constant = 0.95

    with open("scan.inp", "w") as f:
        f.write("$constrain\n")
        f.write("   force constant="+str(force_constant)+"\n")
        f.write("   distance: "+str(args.center_atom[1])+", "
                + str(args.center_atom[0])+", "+str(args.bond1_distance[0])+"\n")
        f.write("$scan\n")
        f.write("   1: "+str(args.bond1_distance[0])+", "+str(args.bond1_distance[1])
                + ", "+str(args.cut_distance[0])+"\n")
        f.write("$opt\n")
        f.write("   maxcycle="+str(args.max_cycle)+"\n")
        f.write("$end")

    subprocess.call("bash 1Atom.sh", shell=True)

    xyzFile = GeometryXYZs("xtbscan.xyz")
    xyzFile.method_read_xyz()
    for idx1 in range(1, len(xyzFile)+1, 1):
        xyzFile.set_filename("xtbscan_single.xyz")

        xyzFile.method_save_xyz([idx1])

        with open("scan.inp", "w") as f:
            # f.write("$fix")
            # f.write("   atoms : "+str(args.center_atom[1])+", "+str(args.center_atom[0])+"\n")
            f.write("$constrain\n")
            f.write("   force constant="+str(force_constant)+"\n")
            f.write("   distance: "+str(args.center_atom[1])+", "
                    + str(args.center_atom[2])+", "+str(args.bond2_distance[0])+"\n")
            f.write("   distance: "+str(args.center_atom[1])+", "
                    + str(args.center_atom[0])+", "+str(abs(args.bond1_distance[0]-args.bond1_distance[1])*(idx1-1)/len(xyzFile)+args.bond1_distance[0])+"\n")
            f.write("$scan\n")
            f.write("   1: "+str(args.bond2_distance[0])+", "+str(args.bond2_distance[1])
                    + ", "+str(args.cut_distance[1])+"\n")
            f.write("$opt\n")
            f.write("   maxcycle="+str(args.max_cycle)+"\n")
            f.write("$end")

        subprocess.call("bash 2Atom.sh", shell=True)
    subprocess.call(
        "grep 'energy' xtbscan2.xyz | awk '{print $2}' > out", shell=True)
    subprocess.call(
        "rm -rf xtbscan_single.xyz xtb_1Atom scan.inp xtb_2Atom", shell=True)
