#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import subprocess
# from icecream import ic
from pathlib import Path

descr = """
________________________________________________________________________________
| For thermo calculation by using GFN-xTB method   
| Usages   : molclus_thermo.py <geometry> [options]
| Input    : -i input file [default traj.xyz]
| [options]
| Method   : --method To set the method gfn0/gfn1/gfn2/gfnff [default gfn2]
| Solvent  : --alpb To set the solvent (for alpb mode and prefered choice)
|             CHCl3/DMSO/H2O(water)  
| Solvent  : --gbsa To set the solvent (for gbsa mode)
|             methanol/CHCl3/DMSO/H2O(water)  
| Charge   : --chrg to set the charge on the molecule [default 0]
| UHF      : --uhf to set the number of unpaired electrons [default 1]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
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
        args = cml()

    inFile = Path(args.file)
    single_xyz_name = Path(".temp.xyz")

    import platform
    _system: str = platform.system()

    # Default to xtb command
    from censo_ext.Tools.utility import program_IsExist
    xtb_cmd: str = ""
    prog = "xtb"
    program_IsExist(prog)
    xtb_cmd += prog

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    from censo_ext.Tools.utility import program_IsExist
    program_IsExist("xtb")

    print(f" Inputted geometry file: {inFile}")
    xtb_cmd += f" {single_xyz_name}"

    print(" Loading basic information from the inputted geometry file ...")
    print(f" There are totally        {len(xyzFile)} geometries in the inputted geometry file\n")  # nopep8
    print(f" Setting method :  {args.method}")
    cmd_solvent = "vacuum"
    if args.alpb:
        cmd_solvent: str = args.alpb
    elif args.gbsa:
        cmd_solvent: str = args.gbsa

    print(f" Setting solvent : {cmd_solvent}")
    print(" Loading setting data ...")
    xtb_cmd += f" --{args.method} --bhess vtight"
    if args.alpb:
        xtb_cmd += f" --alpb {args.alpb}"
    if args.gbsa:
        xtb_cmd += f" --gbsa {args.gbsa}"

    xtb_cmd += f" --chrg {args.chrg} --uhf {args.uhf}"  # nopep8

    print(" All conformer in the inputted geometry file will be processed")
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")

    xtb_cmd += " --enso -I ../xcontrol-inp > thermo.out"
    xcontrol_inp: Path = Path("xcontrol-inp")
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
    for idx1 in range(1, len(xyzFile)+1, 1):
        xyzFile.set_filename(single_xyz_name)
        xyzFile.method_save_xyz([idx1])
        print(f"                          *** Configuration         {idx1}  ****")  # nopep8
        print(f" Loading geometry	 {idx1}  from the inputted geometry file")      # nopep8
        print(" Generating  file...")
        subprocess.call(xtb_cmd, shell=True)
        print(f" Running:  {xtb_cmd}")

        lines: list = open("thermo.out", "r").readlines()
        import re
        for line in lines:
            if re.search(r'contrib\.', line):
                thermo.append(line.split()[3])

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
