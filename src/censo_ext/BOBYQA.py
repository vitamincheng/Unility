#!/usr/bin/env python
from icecream import ic
import argparse
import os
import numpy as np
import numpy.typing as npt
from pathlib import Path
import subprocess
from censo_ext.Tools.anmrfile import AD_BOBYQA, AD_Normal
from censo_ext.Tools.utility import print_descr

descr = """
________________________________________________________________________________
| using BOBYQA method to fit
| Usages   : BOBYQA.py [options]
| [options]
| Dir      : -d the Dir [default .]
| Ref      : -r the actual reference file [default 1r.dat or 1r.npz]
| mf       : -mf magnetic frequency of scan nmr [default 500.0]
| lw       : -lw line width of scan nmr [2.0 for H, 40 for C]
| Limit    : -l limit border(ppm) [defalut 0.20]
| Prog     : -p --prog Use external anmr execute file [default False]
| verbose  : -v --verbose more detail [default False]
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
        "-d",
        "--dir",
        dest="dir",
        action="store",
        required=False,
        default=".",
        help="Provide output_file name [default .]",
    )

    parser.add_argument(
        "-r",
        "--ref",
        dest="ref",
        action="store",
        required=False,
        default="1r.npz",
        help="Provide ref file(dat/npz) name [default 1r.npz]",
    )

    parser.add_argument(
        "-l",
        "--limit",
        dest="limit",
        action="store",
        type=float,
        required=False,
        default=0.20,
        help="Provide limit border (ppm) [0.20]",
    )

    parser.add_argument(
        "-p",
        "--prog",
        dest="prog",
        action="store_true",
        help="Use external anmr execute file [default False]",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose mode [default False]",
    )

    parser.add_argument(
        "-mf",
        "--magnfreq",
        dest="mf",
        action="store",
        type=float,
        required=False,
        default=500.0,
        help="magnetic frequency of scan nmr [default 500.0]",
    )
    parser.add_argument(
        "-lw",
        "-linewidth",
        dest="lw",
        action="store",
        type=float,
        required=False,
        default=2,
        help="line width of scan nmr [default 2.0 for H, 40.0 for C]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


class global_variable():
    Dir: Path
    ref_dat: Path
    limit: float
    prog: bool
    mf: float
    lw: float
    ref: float
    idx_keys: list
    AD_bobyqa: AD_BOBYQA
    AD_normal: AD_Normal

    def __init__(self) -> None:
        self.FileAnmr: Path = Path("anmr.dat")

    def method_update_Dir(self) -> None:
        try:
            self.DirFileAnmr: Path = self.Dir / self.FileAnmr
            self.DirFileRef_dat: Path = self.Dir / self.ref_dat
            self.AD_bobyqa = AD_BOBYQA(g_var.Dir)
            self.AD_normal = AD_Normal(g_var.Dir)
            self.AD_bobyqa.method_load_files()  # all three files are or not exists
            self.AD_normal.method_load_files()  # all three files are or not exists
        except NameError:
            print("  self.Dir is not defined.")
            exit(1)


global g_var
g_var = global_variable()


def rosenbrock(x0: npt.NDArray[np.float64]) -> float:
    """Objective function for BOBYQA optimization of NMR chemical shifts.

    This function updates the chemical shift parameters in ORCA input files,
    runs NMR calculation (either externally via anmr.sh or internally via anmr.py),
    and returns the sum of squared differences between calculated and reference
    NMR spectra.

    Args:
        x0: Array of parameter values to optimize. For single peak optimization,
            this contains one value. For group peak optimization, it contains
            multiple values corresponding to different peaks.

    Returns:
        float: Sum of squared differences between calculated and reference NMR spectra.
    """
    from censo_ext.Tools.datfile import CensoDat
    # SParams_bobyqa: npt.NDArray[np.float64] = np.genfromtxt(
    #    g_var.AD_bobyqa._file_orcaS)
    SParams_bobyqa: npt.NDArray[np.float64] = g_var.AD_bobyqa.SParams  # type: ignore # nopep8

    if len(x0) == 1:
        # single peak
        for idx_key in g_var.idx_keys:
            SParams_bobyqa[idx_key][1] = x0[0]
    else:
        # group peaks
        for idx0, loop in enumerate(g_var.idx_keys):
            for idx_key in loop:  # type: ignore
                SParams_bobyqa[idx_key][1] = x0[idx0]

    g_var.AD_bobyqa.SParams = SParams_bobyqa
    g_var.AD_bobyqa.method_save_files()
    SParams_exec: npt.NDArray[np.float64] = np.delete(
        SParams_bobyqa, 2, axis=1)

    if g_var.prog:
        # print("External program: anmr")
        import sys
        cwd: Path = Path.cwd()
        os.chdir(g_var.Dir)
        template_inp: Path = Path("CONF1/NMR/orcaS.out")
        with open(template_inp, "w") as f:
            sys.stdout = f
            template = """
--------------------------------
CHEMICAL SHIELDING SUMMARY (ppm)
--------------------------------


  Nucleus  Element    Isotropic     Anisotropy
  -------  -------  ------------   ------------
"""
            print(template)
            sys.stdout = sys.__stdout__
        SParams_exec.T[1] = SParams_exec.T[1] + g_var.ref
        SParams_exec.T[0] = SParams_exec.T[0]-1
        file_orcaS_main = Path("CONF1/NMR/orcaS-main.out")
        np.savetxt(file_orcaS_main, SParams_exec,
                   fmt="%7d       H    %10.5f          0")
        subprocess.call(
            f"cat {file_orcaS_main} >> CONF1/NMR/orcaS.out", shell=True)
        subprocess.call(f"rm {file_orcaS_main}", shell=True)

        sys.stdout = open(os.devnull, 'w')
        result = subprocess.call("anmr.sh", shell=True)
        sys.stdout = sys.__stdout__

        if result != 0:
            raise ValueError(" call anmr.sh process have something wrong !!!")  # nopep8

        os.chdir(cwd)
        dat_Sim: CensoDat = CensoDat(file=g_var.DirFileAnmr)

    elif not g_var.prog:
        # print("Internal python: anmr.py")
        g_var.AD_normal.SParams = SParams_exec
        g_var.AD_normal.method_save_files()
        import censo_ext.anmr as anmr
        x: dict = {'out': 'output.npz', "dir": g_var.Dir, "json": None, 'mf': g_var.mf,
                   'lw': g_var.lw, 'ascal': None, 'bscal': None, 'thr': 0.30, 'thrab': 0.020,
                   'tb': 4, 'cutoff': 0.001, 'start': None, 'end': None, "verbose": False,
                   'mss': 10, 'auto': True, 'average': True, 'bobyqa': False}
        import sys
        sys.stdout = open(os.devnull, 'w')
        anmr.main(args=argparse.Namespace(**x))
        sys.stdout = sys.__stdout__

        dat_Sim: CensoDat = CensoDat(file=g_var.Dir / Path(x["out"]))
    else:
        raise ValueError("Something wrong in your argument")

    dat_Sim.method_normalize_dat()
    try:
        dat_diff: CensoDat = dat_Sim - dat_Ref  # type: ignore
    except NameError:
        dat_Ref = CensoDat(file=g_var.DirFileRef_dat)
        dat_Ref.method_normalize_dat()
        dat_diff: CensoDat = dat_Sim - dat_Ref

    return np.sum(np.square(dat_diff.get_Dat()))


def Scan_single_Peak(args) -> None:
    """Perform BOBYQA optimization for individual NMR peaks.

    This function iterates through all unique peak serial numbers in the ORCA
    input file, optimizes each peak individually using the BOBYQA algorithm,
    and prints optimization progress and results.

    Args:
        args: Command-line arguments containing optimization settings.
            Expected attributes include verbose for detailed output.

    Returns:
        None: This function performs optimization but doesn't return a value.
    """
    import pybobyqa

    SParams_bobyqa: npt.NDArray[np.float64] = g_var.AD_bobyqa.SParams.T  # type: ignore # nopep8

    in_sets: set[int] = set(SParams_bobyqa[2].astype(int).tolist())
    in_sets = {x for x in in_sets if x < 1000 and x >= 1}
    print(f"  {in_sets=}")
    print("  ========== Start single_peak ==========")
    for in_set in in_sets:
        print("")
        if args.verbose:
            ic(in_set)
        intp: npt.NDArray[np.int64] = np.argwhere(
            SParams_bobyqa[2] == in_set).flatten()
        SParams: list[float] = list(
            map(float, np.atleast_1d(SParams_bobyqa[1][intp[0]])))
        if args.verbose:
            ic(SParams)
            ic(intp)

        g_var.idx_keys = list(intp)
        x0: npt.NDArray[np.float64] = np.array(SParams)
        bounds = x0 - g_var.limit, x0 + g_var.limit

        print(f"{x0=}")
        soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=bounds,
                              scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
        print(f"{soln.f=} {soln.x=}")
    print("  ========== Finished single_peak ==========")


def Scan_group_Peaks(args) -> None:
    """Perform BOBYQA optimization for groups of coupled NMR peaks.

    This function identifies peak groups (serial numbers >= 1000) in the ORCA
    input file, performs permutation-based optimization of all possible peak
    ordering combinations, and selects the best configuration based on minimum
    objective function value. It then performs a final refinement with tight
    bounds around the best solution.

    Args:
        args: Command-line arguments containing optimization settings.
            Expected attributes include verbose for detailed output.

    Returns:
        None: This function performs optimization but doesn't return a value.
    """
    import pybobyqa
    SParams_bobyqa: npt.NDArray[np.float64] = g_var.AD_bobyqa.SParams.T  # type: ignore # nopep8

    in_sets: set[int] = set(SParams_bobyqa[2].astype(int).tolist())
    in_sets = {x for x in in_sets if x >= 1000}
    if len(in_sets) == 0:
        return

    print("  ========== Start group_peaks ==========")
    print(f"  {in_sets=}")

    # Data structure of idx, Chemical_Shift, idx_atoms
    Data: list[list] = []

    for in_set in in_sets:
        intp: npt.NDArray[np.int64] = np.argwhere(
            SParams_bobyqa[2] == in_set).flatten()
        Data.append(
            [in_set, SParams_bobyqa[1][intp[0]], intp])

    if args.verbose:
        ic(Data)
    nNumbers: int = len(Data)
    from itertools import permutations
    Permutations: list[tuple[int, ...]] = list(permutations(
        [*range(0, nNumbers)], nNumbers))
    solution_f: list = []
    solution_x0: list = []

    for Permutation in Permutations:
        x0: npt.NDArray[np.float64] = np.array([x[1] for x in Data])[
            list(Permutation)]
        print(f"\n  {x0=}")
        g_var.idx_keys = [x[2] for x in Data]
        if args.verbose:
            ic(g_var.idx_keys)
        bounds = x0 - g_var.limit, x0 + g_var.limit
        soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=bounds,
                              scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)
        print(f"{soln.f=} {soln.x=}")
        solution_f.append(soln.f)
        solution_x0.append(soln.x)

    print(f"\n{solution_f=}")
    argsmin: int = min(range(len(solution_f)), key=solution_f.__getitem__)
    list_x0: list[float] = solution_x0[argsmin]

    # After Permutations, the best choice x0 is calculated again and get output.dat or anmr.dat
    print("\n The best choice x0 is calculation again")
    print(f"{list_x0=}")
    x0 = np.array(list_x0)
    limit_tiny: float = 0.0001
    bounds_tiny = x0 - limit_tiny, x0 + limit_tiny
    soln = pybobyqa.solve(rosenbrock, x0, print_progress=True, bounds=bounds_tiny,
                          scaling_within_bounds=True, rhobeg=0.01, rhoend=0.00001)

    print("  ========== Finished group_peaks ==========")


def Create_BOBYQA() -> None:
    """Create initial BOBYQA input file from ORCA output.

    This function reads the original ORCA chemical shift file and creates a
    modified version with an additional column for BOBYQA optimization control.
    The third column determines how each peak will be optimized:
    - Values 1-99: Individual peak optimization
    - Values >= 1000: Group peak optimization (all peaks with same number optimized together)

    Args:
        args: Command-line arguments (not used in this function).

    Returns:
        None: This function creates a file and exits the program.
    """
    SParams: npt.NDArray[np.float64] = np.array(list(g_var.AD_normal.SParams.items()))  # type: ignore # nopep8
    g_var.AD_bobyqa.idx1Atoms = g_var.AD_normal.idx1Atoms
    g_var.AD_bobyqa.JCoups = g_var.AD_normal.JCoups
    g_var.AD_bobyqa.SParams = np.insert(SParams, 2, 0, axis=1)
    g_var.AD_bobyqa.method_save_files()

    descr = """
________________________________________________________________________________
| Create the orcaS-BOBYQA.out file 
| three column :         0 - Do nothing 
|                     1~99 - Use BOBYQA to calcuate and fit each chemical shift 
|                            of each number
|               Above 1000 - Use BOBYQA to calucate and fit each chemical shift 
|                            of all groups to find the peaks in one time
|                            
|                            Each chemical shift of the same number is assigned 
|                            to the same by the first chemical shift
| Run this program again
|_______________________________________________________________________________

"""
    print(descr)
    exit(0)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    # orcaS-BOBYQA.out  for setting
    # Average/NMR/orcaS.out for anmr.py   (internal)
    # CONF1/NMR/orcaS.out for anmr        (external)

    if args == argparse.Namespace():
        args = cml()
    print_descr(descr)

    if args.dir and args.ref:               # default args.dir="." and args.ref="1r.dat"
        g_var.Dir = Path(args.dir)
        g_var.ref_dat = Path(args.ref)
        g_var.method_update_Dir()

    if args.limit:                          # default 0.20 ppm
        g_var.limit = args.limit
    if args.prog:
        g_var.prog = True
    else:
        g_var.prog = False

    if not args.mf or not args.lw:
        print("  Args.mf or args.lw is not set !!!")
        print("  Exit and Close the program !!!")
        exit(1)
    else:
        g_var.mf = args.mf
        g_var.lw = args.lw

    from censo_ext.Tools.anmrfile import Anmr
    inAnmr: Anmr = Anmr(Dir=args.dir, verbose=args.verbose)
    inAnmr.method_read_anmrrc()
    g_var.ref = inAnmr.get_Anmr_Reference_anmrrc()

    if g_var.AD_normal.Exist():
        print(f"  The File {g_var.AD_normal._file_orcaS} is exist")
        if g_var.AD_bobyqa.Exist():
            print(f"  The file {g_var.AD_bobyqa._file_orcaS} is exist")
            if g_var.prog:
                cwd: Path = Path.cwd()
                os.chdir(g_var.Dir)
                print(" Need to build the new CONF* system")
                print(" And copy your CONF* to /Backup/CONF*")
                print(" And create a new CONF1 (copy from /Backup/CONF1)")
                print(" Modify from /Average/NMR/orcaS.out")
                Res = input("Are you Sure to Continue ?? (Y/N)")
                if Res == "Y" or Res == "y":
                    subprocess.call("mkdir backup", shell=True)
                    subprocess.call("mv CONF* backup", shell=True)
                    subprocess.call("cp -r backup/CONF1/ .", shell=True)
                else:
                    print("  Exit and Close the program !!!")
                    exit(0)
                os.chdir(cwd)
            Scan_single_Peak(args)
            Scan_group_Peaks(args)
            if g_var.prog:
                cwd: Path = Path.cwd()
                os.chdir(g_var.Dir)
                subprocess.call("rm -rf CONF1", shell=True)
                subprocess.call("mv backup/CONF* .", shell=True)
                subprocess.call("rmdir backup", shell=True)
                os.chdir(cwd)
        else:
            Create_BOBYQA()
    else:
        raise FileNotFoundError(
            f"{g_var.AD_normal._file_orcaS} is not exist !!!")  # type: ignore # nopep8


if __name__ == "__main__":
    main()

# pytest
#
#
# standard test
# BOBYQA.py -d tests/data/06.EthylAcetate/03.Censo -r 1r.dat
# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) -r 1r_h.dat
#
# Want to use external anmr program (Not Implemented on pytest,because subprocess function)
# Start new one
# Run below anmr.py command and it will create average folder
# anmr.py   -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) --auto
# Run below BOBYQA.py command and it will create ocaS-BOBYQA.out
# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS)
# rewrite orcaS-BOBYQA.out file and run below BOBYQA.py by use external program
# BOBYQA.py -d tests/data/31.Cyclohexanone/03.Censo_For_Hydorgen(revTPSS) -r 1r_h.dat -p
#
