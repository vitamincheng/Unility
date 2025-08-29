#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import subprocess
from pathlib import Path
descr = """
________________________________________________________________________________
|                                          [07.05.2025] vitamin.cheng@gmail.com
| molclus_gxtb.py  
| Usages   : molclus_gxtb.py <geometry> [options]
| [options]
| Input    : -i input file [default traj.xyz]
| Output   : -o output file [default isomers.xyz]
| Opt      : --opt To optimize energy
| Method   : --method To set the method gxtb [default gxtb]
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
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="isomers.xyz",
        help="Provide output_file name [default isomers.xyz]",
    )

    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        required=False,
        default="gxtb",
        help="Method gxtb [default gxtb]",
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
        "--opt",
        dest="opt",
        action="store_true",
        help="Optimize energy",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)

    single_traj_Name = ".single_traj.xyz"
    temp_isomer_Name = ".isomers.xyz"
    xtb_cmd: str = ""
    infile: GeometryXYZs = GeometryXYZs(args.file)
    infile.method_read_xyz()
    from censo_ext.Tools.utility import program_IsExist
    prog = "xtb"
    program_IsExist(prog)
    xtb_cmd += prog

    print(f" Inputted geometry file: {args.file}")
    xtb_cmd += " " + single_traj_Name
    print(" Loading basic information from the inputted geometry file ...")
    print(f" There are totally       {str(len(infile))} geometries in the inputted geometry file\n")  # nopep8
    print(f" Setting method :  {args.method}")
    cmd_solvent = "vacuum"
    if args.alpb:
        cmd_solvent = args.alpb
    elif args.gbsa:
        cmd_solvent = args.gbsa
    print(f" Setting solvent :  {cmd_solvent}")
    print(" Loading setting data ...")
    xtb_cmd += f" --{args.method}"
    if args.alpb:
        xtb_cmd += f" --alpb {args.alpb}"
    if args.gbsa:
        xtb_cmd += f" --gbsa {args.gbsa}"

    xtb_cmd += f" --chrg {str(args.chrg)}  --uhf {str(args.uhf)}"

    print(" All conformer in the inputted geometry file will be processed")
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")

    gxtb_cmd = " --driver \"gxtb -grad -c xtbdriver.xyz\""
    xtb_cmd += gxtb_cmd

    if args.opt:
        # xtb_cmd += " --opt loose"
        xtb_cmd += " --opt"

    xtb_cmd += " > xtb.out"

    for idx in range(1, len(infile)+1, 1):
        # idx_str : str = f"{[idx]:05d}"
        infile.set_filename(Path(single_traj_Name))
        infile.method_save_xyz([idx])
        print(f"                          *** Configuration         {str(idx)}  ****")  # nopep8
        print(f" Loading geometry	 {str(idx)}  from the inputted geometry file")      # nopep8
        print(" Generating  file...")
        subprocess.call(xtb_cmd, shell=True)
        print(f" Running:  {xtb_cmd}")

        get_energy: int | None = None
        if args.opt:
            subprocess.call(
                f"cat xtbopt.xyz >> {temp_isomer_Name}", shell=True)

        else:
            # print("singe point")
            xtb_lines = open("xtb.out", "r").readlines()
            import re
            for idy, y in enumerate(xtb_lines):
                if re.search(r"TOTAL ENERGY", y):
                    get_energy = idy
            if get_energy:
                infile.Sts[idx -
                           1].comment_energy = float(xtb_lines[get_energy].split()[3])

    if args.opt:
        # print("opt")
        outfile: GeometryXYZs = GeometryXYZs(Path(temp_isomer_Name))
        outfile.method_read_xyz()
        outfile.method_comment_new()
        outfile.set_filename(args.out)
        outfile.method_save_xyz([])
    else:
        # print("singe point")
        infile.method_rewrite_comment()
        infile.method_comment_new()
        infile.set_filename(args.out)
        infile.method_save_xyz([])

    subprocess.call(
        "rm -rf charges wbo xtb.out xtbrestart xtbtopo.mol xtbopt* .xtboptok gxtbrestart gradient energy coord xtbdriver.xyz", shell=True)
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(temp_isomer_Name, single_traj_Name)


if __name__ == "__main__":
    main()

#   test
#   python3 molclus_xtb.py -i ../tests/crest_conformers.xyz --alpb CHCl3
