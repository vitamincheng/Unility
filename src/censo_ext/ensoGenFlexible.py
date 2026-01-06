#!/usr/bin/env python
import argparse
import os
import shutil
from pathlib import Path

# from icecream import ic
from censo_ext.Tools.utility import delete_all_files, print_arguments
from censo_ext.Tools.utility import copy_file
from censo_ext.molclus_thermo import thermo_process

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


# global variable
outFile: Path = Path("isomers.xyz")
inFile: Path = Path("traj.xyz")
xtbDir: Path = Path(".xtb")
orcaDir: Path = Path(".orca")
thermoDir: Path = Path(".thermo")


def xtb(_inFile: Path) -> None:
    print(" ========== molclus_xtb.py ==========")
    if not xtbDir.is_dir():
        xtbDir.mkdir()
    copy_file(_inFile, xtbDir / inFile)
    cwd: Path = Path.cwd()
    x: dict = {"file": inFile, "method": "gfn2", "chrg": 0, "uhf": 1,
               "out": outFile, "alpb": "CHCl3", "gbsa": None, "opt": True, "new": False}
    import censo_ext.molclus_xtb as molclus_xtb
    os.chdir(cwd / xtbDir)
    molclus_xtb.main(argparse.Namespace(**x))
    os.chdir(cwd)
    copy_file(xtbDir / outFile, inFile)
    shutil.rmtree(xtbDir, ignore_errors=True)
    print(" ========== End ==========")


def orca() -> None:
    print(" ========== molclus_orca.py ==========")
    if not orcaDir.is_dir():
        orcaDir.mkdir()
    copy_file(inFile, orcaDir / inFile)
    cwd: Path = Path.cwd()

    import censo_ext.molclus_orca as molclus_orca
    x: dict = {"file": inFile, "template": "template.inp", "reserve": False,
               "chrg": 0, "uhf": 1, "out": outFile, "convergence": -1, "new": False}
    os.chdir(cwd/orcaDir)
    molclus_orca.main(argparse.Namespace(**x))
    os.chdir(cwd)
    copy_file(orcaDir / outFile, inFile)
    shutil.rmtree(orcaDir, ignore_errors=True)
    print(" ========== End ==========")


def thermo() -> list[str]:

    print(" ========= molclus_thermo.py ==========")
    if not thermoDir.is_dir():
        thermoDir.mkdir()
    copy_file(inFile, thermoDir / inFile)

    cwd: Path = Path.cwd()
    os.chdir(thermoDir)
    args_x: dict = {"file": inFile, "method": "gfn2",
                    "alpb": "CHCl3", "gbsa": None, "chrg": 0, "uhf": 1}
    thermo: list[str] = thermo_process(argparse.Namespace(**args_x))
    os.chdir(cwd)
    shutil.rmtree(thermoDir, ignore_errors=True)
    print(" ========== End ==========")
    return thermo


def ensoGen(_temp: float, thermo: list[str]) -> None:
    from censo_ext.Tools.xyzfile import GeometryXYZs
    from censo_ext.Tools.anmrfile import Anmr
    print(" ========= ensoGenFlexible ==========")
    xyzFile: GeometryXYZs = GeometryXYZs(inFile)
    xyzFile.method_read_xyz()
    outAnmr: Anmr = Anmr()
    outAnmr.method_create_enso(
        xyzFile.method_ensoGenFlexible(_temp, thermo))
    outAnmr.method_save_enso()
    print(" Saved the anmr_enso.new in your working directory ")
    print(" ========== End ==========")


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()
    print_arguments()

    inFile = Path(args.file)
    # fileName: str = p.name
    # Dir: Path = p.parents[0]

    if args.manual:
        choice_xtb: str = input(
            " Geometry optimization use xtb : Yes or No ").lower().split()[0]
        if choice_xtb == "y" or choice_xtb == "yes":
            xtb(inFile)
        choice_orca: str = input(
            " Geometry optimization use orca : Yes or No ").lower().split()[0]
        if choice_orca == "y" or choice_orca == "yes":
            orca()
        else:
            print(f" Direct use {inFile} as opt xyz file ")
        ensoGen(_temp=args.temp, thermo=thermo())
    else:
        xtb(inFile)
        orca()
        ensoGen(_temp=args.temp, thermo=thermo())
    delete_all_files(inFile, outFile)


if __name__ == "__main__":
    main()
