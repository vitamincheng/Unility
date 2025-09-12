#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
import argparse
import subprocess
from pathlib import Path
descr = """
________________________________________________________________________________
| For GFN-xTB of molecules of xyz file  
| Usages   : molclus_xtb.py <geometry> [options]
| Input    : -i input file [default traj.xyz]
| Output   : -o output file [default isomers.xyz]
| [options]
| Opt      : --opt To optimize energy [default False]
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
        "--opt",
        dest="opt",
        action="store_true",
        help="Optimize energy [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()

    single_traj_Name = Path(".solo.xyz")
    temp_isomer_Name = Path(".isomers.xyz")
    xtb_cmd: str = ""
    inFile = Path(args.file)
    outFile = Path(args.out)
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()

    # Default to xtb command
    from censo_ext.Tools.utility import program_IsExist
    prog = "xtb"
    program_IsExist(prog)
    xtb_cmd += prog

    print(f" Inputted geometry file: {inFile}")
    xtb_cmd += f" {single_traj_Name}"
    print(" Loading basic information from the inputted geometry file ...")
    print(f" There are totally       {len(xyzFile)} geometries in the inputted geometry file\n")  # nopep8
    print(f" Setting method :  {args.method}")
    cmd_solvent = "vacuum"
    if args.alpb and (not args.gbsa):
        cmd_solvent = args.alpb
    elif args.gbsa:
        cmd_solvent = args.gbsa
    print(f" Setting solvent :  {cmd_solvent}")
    print(" Loading setting data ...")
    xtb_cmd += f" --{args.method}"
    if args.alpb and (not args.gbsa):
        xtb_cmd += f" --alpb {args.alpb}"
    elif args.gbsa:
        xtb_cmd += f" --gbsa {args.gbsa}"

    xtb_cmd += f" --chrg {args.chrg}  --uhf {args.uhf}"

    print(" All conformer in the inputted geometry file will be processed")
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")
    if args.opt:
        xtb_cmd += " --opt"

    for idx1 in range(1, len(xyzFile)+1, 1):
        xyzFile.set_filename(single_traj_Name)
        xyzFile.method_save_xyz([idx1])
        print(f"                          *** Configuration         {idx1}  ****")  # nopep8
        print(f" Loading geometry	 {idx1}  from the inputted geometry file")      # nopep8
        print(" Generating  file...")
        subprocess.call(f"{xtb_cmd} > xtb.out", shell=True)
        idx1_str = (f"{idx1:05d}")
        print(f" Running:  {xtb_cmd} > {idx1_str}.out")

        get_energy: int | None = None
        if args.opt:
            subprocess.call(
                f"cat xtbopt.xyz >> {temp_isomer_Name}", shell=True)

        else:
            # print("singe point")
            lines: list[str] = open("xtb.out", "r").readlines()
            import re
            for idy0, y in enumerate(lines):
                if re.search(r"TOTAL ENERGY", y):
                    get_energy = idy0
            if get_energy:
                xyzFile.Sts[idx1 -
                            1].comment_energy = float(lines[get_energy].split()[3])

    if args.opt:
        # print("opt")
        optFile: GeometryXYZs = GeometryXYZs(temp_isomer_Name)
        optFile.method_read_xyz()
        optFile.method_comment_new()
        optFile.set_filename(outFile)
        optFile.method_save_xyz([])
    else:
        # print("singe point")
        xyzFile.method_rewrite_comment()
        xyzFile.method_comment_new()
        xyzFile.set_filename(outFile)
        xyzFile.method_save_xyz([])

    subprocess.call(
        "rm -rf charges wbo xtb.out xtbrestart xtbtopo.mol xtbopt* .xtboptok", shell=True)
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(temp_isomer_Name, single_traj_Name)


if __name__ == "__main__":
    main()

#   test
#   python3 molclus_xtb.py -i ../tests/crest_conformers.xyz --alpb CHCl3
