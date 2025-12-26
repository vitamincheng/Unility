#!/usr/bin/env python
from censo_ext.Tools.xyzfile import GeometryXYZs
import os
import sys
import argparse
import subprocess
from censo_ext.Tools.utility import delete_all_files, IsExist_bool, print_arguments
from pathlib import Path

descr = """
________________________________________________________________________________
| For Orca calculation
| Usages   : molclus.py <geometry> [options]
| Input    : -i input file [default traj.xyz]
| Output   : -o output file [default isomers.xyz]
| [options]
| Template : -t orca template file [default template.inp]
| Reserve  : -r Reserve all files, otherwise will only reserve .gbw .out .xyz three files.
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
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="isomers.xyz",
        help="Provide output_file name [default isomers.xyz]",
    )

    parser.add_argument(
        "-t",
        "--template",
        dest="template",
        action="store",
        required=False,
        default="template.inp",
        help="Provide orca parameter file name default [template.inp]",
    )

    parser.add_argument(
        "-r",
        "--reserve",
        dest="reserve",
        action="store_true",
        help="Reserve all files, otherwise will Only reserve .gbw .out .xyz files [default False]",
    )
    parser.add_argument(
        "-c",
        "--convergence",
        dest="convergence",
        action="store",
        type=int,
        default=0,
        help="Geometry Optimization thresholds : -1/LooseOpt 0/NormalOpt 1/TightOpt 2/VeryTightOpt [default 0]",
    )

    parser.add_argument(
        "--new",
        dest="new",
        action="store_true",
        help="Reordered the serial number of the cluster in xyz file [default False]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile = Path(args.file)
    outFile = Path(args.out)
    solo_xyz: Path = Path("[xyzfile]")
    template_inp: Path = Path(".template.inp")

    # Ensure input file exists
    from censo_ext.Tools.utility import IsExist
    IsExist(inFile)
    template_Exist: bool = IsExist_bool(args.template)

    # Define default template.inp
    if not template_Exist:
        with open(template_inp, "w") as f:
            sys.stdout = f
            print("! r2SCAN-3c miniprint PAL8 CPCM(chloroform) noautostart")
            match args.convergence:
                case -1:
                    print("! LooseOpt")
                case 0:
                    print("! Opt")
                case 1:
                    print("! TightOpt")
                case 2:
                    print("! VeryTightOpt")
                case 100:  # sp: single point
                    pass
            print("* xyzfile 0 1 [xyzfile]")
        sys.stdout = sys.__stdout__

    # Read input file
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    xyzFile_nClusters: list[int] = [x.comment_nClusters for x in xyzFile.Sts]

    # Find orca executable path
    str_env: list[str] = os.environ['PATH'].split(":")
    prog: str = "orca"
    match: list[str] = [x for x in str_env if prog in x]
    if len(match) == 0:
        raise ValueError(f" Need the {prog} Program !!!")

    orca_path = match[0]+f"/{prog}"
#    orca_path="~/Library/orca_6_1_0/orca"
    if template_Exist:
        template_Name: str = str(args.template)[:-4]
    else:
        template_Name: str = str(template_inp)[:-4]

    print(f" Inputted geometry file: {inFile}")
    print(" Loading basic information from the inputted geometry file ...")
    print(
        f" There are totally       {len(xyzFile)} geometries in the inputted geometry file")
    if template_Exist:
        print(f" Setting file : {args.template}")
    else:
        print(" Setting file : use default [r2SCAN-3c / CHCl3] ")
        args.template = str(template_inp)
    print(" Loading setting file ...")
    print(" All conformer in the inputted geometry file will be processed")
    subprocess.call("rm -f isomers.xyz *.tmp", shell=True)
    delete_all_files(solo_xyz)
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")
    templateFile_Exist: bool = False

    for idx1 in range(1, len(xyzFile)+1, 1):
        idx1_str = (f"{idx1:05d}")
        xyzFile.set_filename(solo_xyz)
        xyzFile.method_save_xyz([idx1])

        print(f"                          "
              f"*** Configuration        {idx1} ****")
        print(f" Loading geometry	{idx1} from the inputted geometry file")
        print(" Generating  file...")

        # Run orca
        orca_cmd: str = f"{orca_path} {args.template} > {template_Name}.out"
        subprocess.call(orca_cmd, shell=True)
        print(f" Running: {orca_path} {args.template} > {idx1_str}.out")

        orca_lines: list[str] = open(template_Name + ".out", "r").readlines()
        get_energy: int | None = None
        for idy0, line in enumerate(orca_lines):
            if r"FINAL SINGLE POINT ENERGY" in line:
                get_energy = idy0

        from os.path import exists
        templateFile_Exist = exists(f"{template_Name}.xyz")
        if templateFile_Exist:
            templateLines: list[str] = open(f"{template_Name}.xyz", "r").readlines()  # nopep8
            for idy0, y in enumerate(templateLines):
                if rf"Coordinates from ORCA-job {template_Name}" in y and get_energy:
                    # get_comment_template = idy
                    templateLines[idy0] = str(
                        orca_lines[get_energy].split()[4] + "\n")
            open(f"{template_Name}.xyz", "w").writelines(templateLines)

            subprocess.call(
                f"cat {template_Name}.xyz >> {str(outFile)}", shell=True)
            subprocess.call(
                f"mv -f {template_Name}.xyz {idx1_str}.xyz", shell=True)
        else:
            if get_energy:
                xyzFile.Sts[idx1 - 1]._comment_energy = float(orca_lines[get_energy].split()[4])  # nopep8

        subprocess.call(f"mv -f {template_Name}.out {idx1_str}.out", shell=True)  # nopep8
        subprocess.call(f"mv -f {template_Name}.gbw {idx1_str}.gbw", shell=True)  # nopep8

    if templateFile_Exist:  # template File is Exists
        optFile: GeometryXYZs = GeometryXYZs(outFile)
        optFile.method_read_xyz()

        for a, b in zip(optFile.Sts, xyzFile_nClusters):
            a.comment_nClusters = b
        optFile.method_rewrite_comment()
        if args.new:
            optFile.method_comment_new()

        optFile.method_save_xyz([])
        print(f" Saved to  {outFile} \n All is done !!!")

    else:
        if args.new:
            xyzFile.method_comment_new()

        xyzFile.set_filename(outFile)
        xyzFile.method_save_xyz([])
        print(f" Saved to  {outFile} \n All is done !!!")

    if not args.reserve:
        subprocess.call(
            f"rm -rf {template_Name}.cpcm {template_Name}.densities {template_Name}.engrad {template_Name}.out {template_Name}_property.txt {template_Name}_trj.xyz {template_Name}.opt", shell=True)
        subprocess.call(
            f"rm -rf {template_Name}.cpcm_corr {template_Name}.densitiesinfo {template_Name}.property.txt {template_Name}.bibtex", shell=True)
        delete_all_files(solo_xyz, template_inp)


if __name__ == "__main__":
    main()

#   test
#   python3 molclus_orca.py -i ../tests/crest_conformers1.xyz (use default)
