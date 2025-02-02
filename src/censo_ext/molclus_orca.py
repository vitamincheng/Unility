#!/usr/bin/env python3
from censo_ext.Tools.xyzfile import GeometryXYZs
import os
import sys
import argparse
import subprocess
from censo_ext.Tools.utility import delete_all_files, is_exist_return_bool
from pathlib import Path

descr = """
________________________________________________________________________________
|                                          [08.27.2024] vitamin.cheng@gmail.com
| Simple crest of multixyz file for orca
| Usages   : molclus.py <geometry> [options]
| [options]
| Input    : -i input file [default traj.xyz]
| Output   : -o output file [default isomers.xyz]
| Template : -t orca template file [default template.inp]
| Remove   : -r only reserve .gbw .out .xyz three files, others will be removed  
| Packages : Tools 
| Module   : xyzfile.py / unility.py
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
        "--remove",
        dest="remove",
        action="store_true",
        required=False,
        default=True,
        help="Only reserve .gbw .out .xyz files, others will be reomved [default True]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)

    from censo_ext.Tools.utility import is_exist
    is_exist(args.file)
    IsExist_template: bool = is_exist_return_bool(args.template)

    template_inp: str = ".template.inp"
    if IsExist_template == False:
        original_stdout = sys.stdout
        with open(template_inp, "w") as f:
            sys.stdout = f
            print("! r2SCAN-3c opt miniprint PAL8 CPCM(chloroform) noautostart")
            print("%maxcore 6000")
            print("* xyzfile 0 1 [xyzfile]")
        sys.stdout = original_stdout

    infile: GeometryXYZs = GeometryXYZs(args.file)
    infile.method_read_xyz()
    singleXYZ: Path = Path("[xyzfile]")

    str_list = os.environ['PATH'].split(":")
    match = [x for x in str_list if "orca" in x]
    print(match)
    if len(match) == 0:
        print(" Need the orca program !!!!")
        print(" Exit the program ")
        exit(0)
    orca_path = match[0]+"/orca"

#    orca_path="~/orca_5_0_4_linux_x86-64_shared_openmpi411/orca"
#    orca_path="~/orca_6_0_0_linux_x86-64_avx2_shared_openmpi416/orca"
    if IsExist_template == True:
        templateName: str = args.template[:-4]
    else:
        templateName: str = template_inp[:-4]

    print(" Inputted geometry file: " + args.file)
    print(" Loading basic information from the inputted geometry file ...")
    print(" There are totally       " + str(len(infile)) +
          " geometries in the inputted geometry file\n")
    if IsExist_template == True:
        print(" Setting file : " + args.template)
    else:
        print(" Setting file : use default [r2SCAN-3c / CHCl3] ")
        args.template = template_inp
    print(" Loading setting file ...")
    print(" All conformer in the inputted geometry file will be processed")
    subprocess.call("rm -f isomers.xyz *.tmp", shell=True)
    delete_all_files(singleXYZ)
    print(" Cleaning old input and temporary files ...")
    print(" Running: rm isomers.xyz *.tmp")
    templateFileIsExists: bool = False

    for idx1 in range(1, len(infile)+1, 1):
        idx1_str = ("{:05d}".format(idx1))
        infile.set_filename(singleXYZ)
        infile.method_save_xyz([idx1])

        print("                          *** Configuration        "+str(idx1)+" ****")
        print(" Loading geometry	"+str(idx1) +
              " from the inputted geometry file")
        print(" Generating  file...")

        orca_cmd: str = orca_path + " "+args.template+" > " + templateName + ".out "
        subprocess.call(orca_cmd, shell=True)
        print(" Running: " + orca_cmd)

        orca_lines: list[str] = open(templateName + ".out", "r").readlines()

        import re
        get_energy: int | None = None
        for idy, y in enumerate(orca_lines):
            if re.search(r"FINAL SINGLE POINT ENERGY", y):
                get_energy = idy

        from os.path import exists
        templateFileIsExists = exists(templateName + ".xyz")
        if templateFileIsExists:
            templateLines: list = open(
                templateName + ".xyz", "r").readlines()
            for idy, y in enumerate(templateLines):
                if re.search(r"Coordinates from ORCA-job "+templateName, y) and get_energy:
                    # get_comment_template = idy
                    templateLines[idy] = str(
                        orca_lines[get_energy].split()[4] + "\n")
            open(templateName + ".xyz", "w").writelines(templateLines)

            subprocess.call("cat " + templateName +
                            ".xyz >> " + args.out, shell=True)
            subprocess.call("mv -f " + templateName +
                            ".xyz " + idx1_str + ".xyz", shell=True)
        else:
            if get_energy:
                infile.Sts[idx1 -
                           1].comment_energy = float(orca_lines[get_energy].split()[4])

        subprocess.call("mv -f " + templateName + ".out " +
                        idx1_str+".out", shell=True)
        subprocess.call("mv -f " + templateName + ".gbw " +
                        idx1_str+".gbw", shell=True)

    if templateFileIsExists:
        outfile: GeometryXYZs = GeometryXYZs(args.out)
        outfile.method_read_xyz()
        outfile.method_comment_new()
        outfile.method_save_xyz([])
        print(f" Saved to  {args.out} \n All is done !!!")
    else:
        infile.method_rewrite_comment()
        infile.method_comment_new()
        infile.set_filename(args.out)
        infile.method_save_xyz([])
        print(f" Saved to  {args.out} \n All is done !!!")

    if args.remove == True:
        subprocess.call("rm -rf "+templateName+".cpcm "+templateName+".densities "+templateName+".engrad "+templateName +
                        ".out "+templateName+"_property.txt "+templateName+"_trj.xyz "+templateName+".opt ", shell=True)
        subprocess.call("rm -rf "+templateName+".cpcm_corr "+templateName+".densitiesinfo "+templateName+".property.txt "+templateName +
                        ".bibtex ", shell=True)
        delete_all_files(singleXYZ, template_inp)


if __name__ == "__main__":
    main()

#   test
#   python3 molclus_orca.py -i ../tests/crest_conformers1.xyz (use default)
