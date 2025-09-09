#!/usr/bin/env python
import argparse
import os
import shutil
from pathlib import Path
# import numpy as np
# from icecream import ic
from censo_ext.Tools.utility import delete_all_files
from censo_ext.Tools.utility import copy_file
from sys import argv as sysargv
import filecmp

descr = """
________________________________________________________________________________
| ensoGenFlexible.py
| Usages   : ensoGenFlexible.py <geometry> [options]
| Input    : -i xyz file [default traj.xyz]
| Output   : -o output anmr_enso file [default anmr_enso.flexible]
| [options]
| Manual   : -m Manually assign the splitting position [default False]
| Temp     : -t the temperature of the environment [default 298.15 K]
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
        help="Provide one input xyz file [default traj.xyz]",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="output.xyz",
        help="Provide one anmr_enso file [default anmr_enso.flexible]",
    )

    parser.add_argument(
        "-m",
        "--manual",
        dest="manual",
        action="store_true",
        help="Assign the splitting position of Atoms [static Atoms, rotation Atoms] [default False]",
    )
    parser.add_argument(
        "-t",
        "--temp",
        dest="temp",
        action="store",
        required=False,
        default=298.15,
        type=float,
        help="Degrees of Temperature [defalut 298.15 K]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


######
outFile: Path = Path("isomers.xyz")
inFile: Path = Path("traj.xyz")


def xtb(args):
    print(" ========== molclus_xtb.py ==========")
    xtb_folder: Path = Path(".xtb")
    if not xtb_folder.is_dir():
        xtb_folder.mkdir()
    copy_file(args.file, xtb_folder / inFile)
    cwd: Path = Path(os.getcwd())
    x: dict = {"file": inFile, "method": "gfn2", "chrg": 0, "uhf": 1,
               "out": outFile, "alpb": "CHCl3", "gbsa": None, "opt": True}
    import censo_ext.molclus_xtb as molclus_xtb
    os.chdir(cwd/xtb_folder)
    molclus_xtb.main(argparse.Namespace(**x))
    os.chdir(cwd)
    copy_file(xtb_folder / outFile, inFile)
    shutil.rmtree(xtb_folder, ignore_errors=True)
    print(" ========== End ==========")


def orca(args, Dir, FileName):
    print(" ========== molclus_orca.py ==========")
    orca_folder: Path = Path(".orca")
    if not orca_folder.is_dir():
        orca_folder.mkdir()
    copy_file(Path(inFile), orca_folder / inFile)
    cwd: Path = Path(os.getcwd())

    import censo_ext.molclus_orca as molclus_orca
    x: dict = {"file": inFile, "template": "template.inp", "reserve": False,
               "chrg": 0, "uhf": 1, "out": outFile}
    os.chdir(cwd/orca_folder)
    molclus_orca.main(argparse.Namespace(**x))
    os.chdir(cwd)
    copy_file(orca_folder / outFile, inFile)
    shutil.rmtree(orca_folder, ignore_errors=True)
    print(" ========== End ==========")


def thermo(args) -> list:
    import censo_ext.molclus_thermo as molclus_thermo
    print(" ========= molclus_thermo.py ==========")
    thermo_folder: Path = Path(".thermo")
    if not thermo_folder.is_dir():
        thermo_folder.mkdir()
    copy_file(inFile, thermo_folder / inFile)

    cwd: Path = Path(os.getcwd())
    os.chdir(thermo_folder)
    args_x: dict = {"file": inFile, "method": "gfn2",
                    "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1}
    thermo: list = molclus_thermo.main(argparse.Namespace(**args_x))
    os.chdir(cwd)
    shutil.rmtree(thermo_folder, ignore_errors=True)
    print(" ========== End ==========")
    return thermo


def ensoGen(args, thermo_list) -> None:
    from censo_ext.Tools.xyzfile import GeometryXYZs
    from censo_ext.Tools.anmrfile import Anmr
    print(" ========= ensoGenFlexible ==========")
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    outAnmr: Anmr = Anmr()
    outAnmr.method_create_enso(
        xyzFile.method_ensoGenFlexible(args, thermo_list))
    outAnmr.method_save_enso()
    print(" Saved the anmr_enso.new in your working directory ")
    print(" ========== End ==========")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print(f"    provided arguments: {" ".join(sysargv)}")

    p = Path(args.file)
    fileName: str = p.name
    Dir: Path = p.parents[0]

    if args.manual:
        choice_xtb: str = input(
            " Geometry optimization use xtb : Yes or No ").lower().split()[0]
        if choice_xtb == "y" or choice_xtb == "yes":
            xtb(args)
            # molclus_xtb.py
        choice_orca: str = input(
            " Geometry optimization use orca : Yes or No ").lower().split()[0]
        if choice_orca == "y" or choice_orca == "yes":
            orca(args, Dir, fileName)
        else:
            print(f" Direct use {args.file} as opt xyz file ")
        thermo_list: list = thermo(args)
        ensoGen(args, thermo_list)
    else:
        xtb(args)
        orca(args, Dir, fileName)
        thermo_list = thermo(args)
        ensoGen(args, thermo_list)
    delete_all_files(inFile, outFile)


if __name__ == "__main__":
    main()
