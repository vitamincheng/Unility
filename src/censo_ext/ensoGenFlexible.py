#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import numpy as np
from icecream import ic
from Tools.utility import Delete_All_files
from Tools.xyzfile import ClassGeometryXYZs
from sys import argv as sysargv

descr = """
________________________________________________________________________________
|                                          [11.08.2024] vitamin.cheng@gmail.com
| ensoGenFlexible.py
| Usages   : ensoGenFlexible.py <geometry> [options]
| [options]
| Input    : -i xyz file [default traj.xyz]
| Output   : -o output anmr_enso file [default anmr_enso.flexible]
| Manual   : -m Manually assign the splitting position [default False]
| Temp     : -t the temperature of the environment [default 298.15 K]
| Packages : Tools
| Module   : xyzfile.py / anmrfile.py /molclus_thermo.py / molclus_orca.py
|            molclus_thermo.py
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
        help="Provide one anmr_enso file [default anmr_enso.flexible]",
    )

    parser.add_argument(
        "-m",
        "--manual",
        dest="manual",
        action="store_true",
        required=False,
        default=False,
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


#####
# G
out_file: Path = Path(".isomers.xyz")
in_file: Path = Path(".traj.xyz")


def xtb(args):
    print(" ========== molclus_xtb.py ==========")
    from Tools.utility import Copy_file
    import os
    temp_folder: Path = Path(".xtb")
    # if not os.path.isdir(temp_folder):
    if not temp_folder.is_dir():
        temp_folder.mkdir()
    Copy_file(args.file, temp_folder / in_file)

    working_Dir = os.getcwd()
    temp_folder.mkdir()
    x: dict = {"file": in_file, "method": "gfn2",
               "chrg": 0, "uhf": 1, "out": out_file, "alpb": "CHCl3", "gbsa": None, "opt": True}
    import molclus_xtb
    molclus_xtb.main(argparse.Namespace(**x))
    os.chdir(working_Dir)
    Copy_file(temp_folder / out_file, in_file)
    Copy_file(temp_folder / out_file, Path("isomers.xyz"))
    print(" Saved the isomers.xyz in your working directory ")
    import shutil
    shutil.rmtree(temp_folder, ignore_errors=True)
    print(" ========== End ==========")


def orca(args, Dir, FileName):
    print(" ========== molclus_orca.py ==========")
    from Tools.utility import Copy_file
    import os
    temp_folder: Path = Path(".orca")
    if not temp_folder.is_dir():
        temp_folder.mkdir()
    Copy_file(Path(in_file), temp_folder / in_file)

    working_Dir = os.getcwd()
    temp_folder.mkdir()
    import molclus_orca

    x: dict = {"file": in_file, "template": "template.inp", "remove": True,
               "chrg": 0, "uhf": 1, "out": out_file}
    molclus_orca.main(argparse.Namespace(**x))
    os.chdir(working_Dir)
    Copy_file(temp_folder / out_file, in_file)
    Copy_file(temp_folder / out_file, Path("isomers.xyz"))
    print(" Saved the isomers.xyz in your working directory ")
    import shutil
    shutil.rmtree(temp_folder, ignore_errors=True)
    print(" ========== End ==========")


def thermo(args) -> list:
    import molclus_thermo
    print(" ========= molclus_thermo.py ==========")
    from Tools.utility import Copy_file
    import os
    thermo_folder: Path = Path(".thermo")
    if not thermo_folder.is_dir():
        thermo_folder.mkdir()
    Copy_file(in_file, thermo_folder / in_file)

    working_Dir = os.getcwd()
    os.chdir(thermo_folder)
    args_x: dict = {"file": in_file,
                    "method": "gfn2", "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1}
    thermo: list = molclus_thermo.main(argparse.Namespace(**args_x))
    os.chdir(working_Dir)
    import shutil
    shutil.rmtree(thermo_folder, ignore_errors=True)
    print(" ========== End ==========")
    return thermo


def ensoGen(args, thermo_list) -> None:
    from Tools.xyzfile import ClassGeometryXYZs
    from Tools.anmrfile import ClassAnmr
    print(" ========= ensoGenFlexible ==========")
    xyzfile: ClassGeometryXYZs = ClassGeometryXYZs(args.file)
    xyzfile.method_read_xyz()
    # outAnmr: ClassAnmr = ClassAnmr(Path(args.file).parents[0])
    outAnmr: ClassAnmr = ClassAnmr()
    outAnmr.method_create_enso(
        xyzfile.method_ensoGenFlexible(args, thermo_list))
    outAnmr.method_save_enso()
    print(" Saved the anmr_enso.new in your working directory ")
    print(" ========== End ==========")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
    print(descr)  # Program description
    print("    provided arguments: {}".format(" ".join(sysargv)))

    WorkingDir: str = os.getcwd()
    # from Tools.utility import IsExistsDirFileName
    args.file
    p = Path(args.file)
    fileName = p.name
    Dir = p.parents[0]
    # Dir, Name = IsExistsDirFileName(Path(args.dir) / Path(args.file))
    # args.file, args.dir = Name, str(Dir)

    if args.manual == True:
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
    Delete_All_files(in_file, out_file)


if __name__ == "__main__":
    main()
