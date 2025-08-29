#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import subprocess
# from icecream import ic
from pathlib import Path

descr = """
________________________________________________________________________________
|                                          [11.07.2024] vitamin.cheng@gmail.com
| molclus_thermo.py  
| Usages   : molclus_thermo.py <geometry> [options]
| [options]
| Input    : -i input file [default traj.xyz]
| Method   : --method To set the method gfn0/gfn1/gfn2/gfnff [default gfn2]
| Solvent  : --alpb To set the solvent (for alpb mode and prefered choice)
|             CHCl3/DMSO/H2O(water)  
| Solvent  : --gbsa To set the solvent (for gbsa mode)
|             methanol/CHCl3/DMSO/H2O(water)  
| Charge   : --chrg to set the charge on the molecule [default 0]
| UHF      : --uhf to set the number of unpaired electrons [default 1]
| Packages : Tools 
| Module   : xyzfile.py / unility.py 
|______________________________________________________________________________
"""


def cml(descr) -> argparse.Namespace:
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
        help="Provide input_file name [default traj.xyz]",
    )

    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        required=False,
        default="gfn2",
        help="Method : gfn0/gfn1/gfn2/gfnff [default gfn2]",
    )

    parser.add_argument(
        "--alpb",
        dest="alpb",
        action="store",
        required=False,
        help="Provide used the solvnet CHCl3/DMSO/H2O(water) ",
    )

    parser.add_argument(
        "--gbsa",
        dest="gbsa",
        action="store",
        required=False,
        help="Provide used the solvnet methanol/CHCl3/DMSO/H2O(water) ",
    )

    parser.add_argument(
        "--chrg",
        dest="chrg",
        action="store",
        type=int,
        required=False,
        default=0,
        help="to set the charge on the molecule [default 0]",
    )

    parser.add_argument(
        "--uhf",
        dest="uhf",
        action="store",
        type=int,
        required=False,
        default=1,
        help="to set the number of unpaired electrons [default 1]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> list[str]:

    if args == argparse.Namespace():
        args = cml(descr)

    single_xyz_name = ".temp.xyz"

    import platform
    _system: str = platform.system()

    # Default to xtb command
    xtb_cmd: str = "xtb"
    if _system == "Darwin":
        xtb_cmd = "/opt/homebrew/bin/xtb"

    infile: GeometryXYZs = GeometryXYZs(args.file)
    infile.method_read_xyz()

    from censo_ext.Tools.utility import program_IsExist
    program_IsExist("xtb")

    print(" Inputted geometry file: "+args.file)
    xtb_cmd = xtb_cmd + " " + single_xyz_name

    print(" Loading basic information from the inputted geometry file ...")
    print(f" There are totally        {str(len(infile))} geometries in the inputted geometry file\n")  # nopep8
    print(f" Setting method :  {args.method}")
    cmd_solvent = "vacuum"
    if args.alpb:
        cmd_solvent: str = args.alpb
    elif args.gbsa:
        cmd_solvent: str = args.gbsa

    print(" Setting solvent : " + cmd_solvent)
    print(" Loading setting data ...")
    xtb_cmd = xtb_cmd + " --" + args.method + " --bhess vtight"
    if args.alpb:
        xtb_cmd = xtb_cmd + " --alpb " + args.alpb
    if args.gbsa:
        xtb_cmd = xtb_cmd + " --gbsa " + args.gbsa

    xtb_cmd = xtb_cmd + " --chrg " + str(args.chrg) + " --uhf " + str(args.uhf)  # nopep8

    print(" All conformer in the inputted geometry file will be processed")
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")

    xtb_cmd = xtb_cmd + " --enso -I ../xcontrol-inp > thermo.out"
    xcontrol_inp: str = "xcontrol-inp"
    import sys
    with open(xcontrol_inp, "w") as f:
        sys.stdout = f
        print("$thermo")
        print("    temp=298.15")
        print("    sthr=50.0")
        print("    imagthr=-100")
        print("$symmetry")
        print("     maxat=1000")
        print("$gbsa")
        print("  gbsagrid=tight")
        print("$end")
    sys.stdout = sys.__stdout__

    thermo: list = []
    for idx in range(1, len(infile)+1, 1):
        infile.set_filename(Path(single_xyz_name))
        infile.method_save_xyz([idx])
        print(f"                          *** Configuration         {str(idx)}  ****")  # nopep8
        print(f" Loading geometry	 {str(idx)}  from the inputted geometry file")      # nopep8
        print(" Generating  file...")
        subprocess.call(xtb_cmd, shell=True)
        print(f" Running:  {xtb_cmd}")

        thermo_lines: list = open("thermo.out", "r").readlines()
        import re
        for idy, y in enumerate(thermo_lines):
            if re.search(r'contrib\.', y):
                thermo.append(y.split()[3])

    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(single_xyz_name, xcontrol_inp)

    return thermo


if __name__ == "__main__":
    main()


#   test
#   python3 molclus_thermo.py -i ../tests/crest_conformers1.xyz --alpb CHCl3 (xcontrol-inp is default, see 148 lines)
#   xtb ../tests/crest_conformers1.xyz --gfn2 --bhess vtight --alpb CHCl3 --chrg 0 --uhf 1 --enso -I xcontrol-inp
#
#
#   xcontrol-inp [fileName]
#   $thermo
#       temp=298.15
#       sthr=50.0
#       imagthr=-100
#   $symmetry
#        maxat=1000
#   $gbsa
#     gbsagrid=tight
#   $end
