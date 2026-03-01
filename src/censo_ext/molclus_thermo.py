#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
from censo_ext.Tools.utility import print_arguments
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

    parser.add_argument(
        "--enso",
        dest="enso",
        action="store_true",
        help="Create the anmr_enso files [default False]",
    )

    parser.add_argument(
        "-t",
        "--temp",
        dest="temp",
        action="store",
        type=float,
        default=298.15,
        help="set the temperature degree K [default 298.15 K]",
    )

    return parser.parse_args()


def thermo_process(args) -> list[str]:
    inFile = Path(args.file)
    single_xyz_name = Path(".temp.xyz")

    import platform
    _system: str = platform.system()

    # Default to xtb command
    from censo_ext.Tools.utility import prog_IsExist
    xtb_cmd: str = ""
    prog = "xtb"
    prog_IsExist(prog)
    xtb_cmd += prog

    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    from censo_ext.Tools.utility import prog_IsExist
    prog_IsExist("xtb")

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
    with open(xcontrol_inp, "w") as f:
        print("$thermo", file=f)
        print("    temp=298.15", file=f)
        print("    sthr=50.0", file=f)
        print("    imagthr=-100", file=f)
        print("$symmetry", file=f)
        print("     maxat=1000", file=f)
        print("$gbsa", file=f)
        print("  gbsagrid=tight", file=f)
        print("$end", file=f)

    thermo: list[str] = []
    entropy: list[str] = []
    for idx1 in range(1, len(xyzFile)+1, 1):
        xyzFile.set_filename(single_xyz_name)
        xyzFile.method_save_xyz([idx1])
        print(f"                          *** Configuration         {idx1}  ****")  # nopep8
        print(f" Loading geometry	 {idx1}  from the inputted geometry file")      # nopep8
        print(" Generating  file...")
        subprocess.call(xtb_cmd, shell=True)
        print(f" Running:  {xtb_cmd}")

        lines: list = open("thermo.out", "r").readlines()
        for line in lines:
            if r'G(RRHO) contrib.' in line:
                thermo.append(line.split()[3])
            if r'TOT   ' in line:
                entropy.append(line.split()[3])

    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(single_xyz_name, xcontrol_inp)
    delete_all_files("charges", "g98.out", "hessian", "thermo.out")
    delete_all_files("vibspectrum", "wbo", "xtb_enso.json")
    delete_all_files("xtbopt.log", "xtbopt.xyz",
                     "xtbrestart", "xtbtopo.mol", "xtbhess.xyz")
    print(entropy)
    import numpy as np
    print("Average Entropy: ", np.average(
        np.array([float(x) for x in entropy])))
    print(thermo)

    return thermo


def enso_generate(args: argparse.Namespace, thermo: list[str], temp: float) -> None:

    from censo_ext.Tools.xyzfile import GeometryXYZs
    from censo_ext.Tools.anmrfile import Anmr
    print("  ===== create anmr_enso =====")
    xyzFile: GeometryXYZs = GeometryXYZs(args.file)
    xyzFile.method_read_xyz()
    outAnmr: Anmr = Anmr()
    outAnmr.method_create_enso(
        xyzFile.method_ensoGenFlexible(temp, thermo))
    outAnmr.method_save_enso()
    print(" Saved the anmr_enso.new in your working directory ")
    print("  ===== End =====")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    thermo: list[str] = thermo_process(args)

    if args.enso and args.temp:
        enso_generate(args=args, thermo=thermo, temp=args.temp)


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
