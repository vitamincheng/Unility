#!/usr/bin/env python
import argparse
import numpy as np
import numpy.typing as npt
import math
# from icecream import ic
from sys import argv as sysargv
from pathlib import Path

descr = """
________________________________________________________________________________
| Usage   : ensoAnalysis.py <enso file> [Options]
| Input   : -i input file [default anmr_enso]
| [Options]
| New     : -n -new copy your inupt file to backup file for new project
| Temp    : -t temperature K [default: 298.15 K]
| Switch  : -s switch ONOFF in anmr_enso file
| Weight  : -w For weights(precent) for every CONFS
| Verbose : -v Verbose mode and show the detail of calculation [default False]
| Output  : -o output file [default average_enso]
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
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="average_enso",
        help="Provide name of the output file without file ending. [average_enso]",
    )
    parser.add_argument(
        "-n",
        "--new",
        dest="new",
        action="store_true",
        help="New file and cpoy input_file to input_file.backup for reference Energy [default False]",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="anmr_enso",
        help="Provide input_file name [default anmr_enso]",
    )
    parser.add_argument(
        "-w",
        "--weights",
        dest="weights",
        action="store_true",
        help="Calculated and Saved the percentage of weight of CONFS to original file [default False]",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Verbose mode and show the detail of Calculation [default False]",
    )

    parser.add_argument(
        "-t",
        "--temperature",
        dest="temp",
        action="store",
        required=False,
        #        default=298.15,
        type=float,
        help="Degrees of Temperature [defalut 298.15 K]",
    )
    parser.add_argument(
        "-s",
        "--switch",
        dest="switch",
        action="store",
        required=False,
        nargs="+",
        type=int,
        help="ONOFF in anmr_enso",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def Boltzmann_enso(np_enso: npt.NDArray, TEMP) -> npt.NDArray:
    import numpy.lib.recfunctions as rfn
    from censo_ext.Tools.Parameter import Eh, FACTOR
    # dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'),
    #       ('Energy', '<f8'), ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8')]

    # Column 8 is Total Gibbs Free Energy (Eh) = Energy + mRRHO
    Total: npt.NDArray[np.float64] = np.array(
        (np_enso['Energy']+np_enso['mRRHO']), dtype=[('Total', 'f8')])
    np_enso = rfn.merge_arrays((np_enso, Total), flatten=True)

    # Gibbs_min is lowest energy of Gibbs Free Energy
    Gibbs_min: np.float64 = np_enso['Total'].min()

    # Column 9 is delta Gibbs Free Energy (kcal/mol)
    Gibbs: npt.NDArray[np.float64] = np.array(
        (np_enso['Total']-Gibbs_min)*Eh, dtype=[('Gibbs', 'f8')])
    np_enso = rfn.merge_arrays((np_enso, Gibbs), flatten=True)

    # Column 1o is Qi (each CONFS)
    Qi: npt.NDArray[np.float64] = np.array(
        np.exp(-np_enso['Gibbs']/(TEMP*FACTOR)), dtype=[('Qi', 'f8')])
    np_enso = rfn.merge_arrays((np_enso, Qi), flatten=True)

    # Qall is sum of Qi
    Qall: np.float64 = np.sum(np_enso['Qi'])

    # Column 11 is percentage of each CONFS
    NEW_BW: npt.NDArray[np.float64] = np.array(
        np_enso['Qi']/Qall, dtype=[('NEW_BW', 'f8')])

    return rfn.merge_arrays((np_enso, NEW_BW), flatten=True)


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml(descr)
        print(descr)  # Program description
        print(f"    provided arguments: {" ".join(sysargv)}")

    from censo_ext.Tools.utility import IsExist_return_bool
    from censo_ext.Tools.Parameter import Eh, Rcal

    inFile: Path = Path(args.file)
    backupFile: Path = Path(str(args.file) + ".backup")

    fileExists: bool = IsExist_return_bool(inFile)
    backupfileExists: bool = IsExist_return_bool(backupFile)

    print("")
    print(f" Reading the input file  : {inFile}")
    print(f" Reading the backup file : {backupFile}")

    if fileExists:
        if backupfileExists:
            pass
        else:

            print(f" The backup file is not exist. {backupFile}")
            print(f" ONOFF args.new : {args.new}")
            if args.new:
                print(
                    f" Copy {inFile} to {backupFile} for original Energy and for reference")
                from censo_ext.Tools.utility import copy_file
                copy_file(inFile, backupFile)
                print("  Run this program again ")
                print("  Exit and Close the program !!!")
                exit(0)

            else:
                print("  Something is Wrong !!!")
                print("  Input file is exists but backup file is not exist. ")
                print(
                    "  IF you want to Create a New project, please Add -n or --new in arugment.")
                print("  Exit and Close the program !!!")
                exit(0)

    else:
        raise FileNotFoundError(
            f"{inFile} was not found or is a directory")

    if (not fileExists) or (not backupfileExists):
        print(f"    {inFile} or {backupFile} , the file is not exist ...")
        print("  Exit and Close the program !!!")
        exit(0)

    np.set_printoptions(precision=4, linewidth=75, suppress=True)

    # anmr_enso    main cal.
    # backup_eneo  as reference
    # a = np.rec.array(np.genfromtxt("limit.out", names=True,
    #                 dtype=['i8', 'i8', 'f8', 'f8']))
    anmr_enso: npt.NDArray = np.genfromtxt(args.file, names=True, dtype=[
        ('i8'), ('i8'), ('i8'), ('f8'), ('f8'), ('f8'), ('f8'), ('f8')])
    backup_enso: npt.NDArray = np.genfromtxt(backupFile, names=True, dtype=[
        ('i8'), ('i8'), ('i8'), ('f8'), ('f8'), ('f8'), ('f8'), ('f8')])

    # dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'),
    #       ('Energy', '<f8'), ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8')]

    if Path("average_enso").exists():
        Path("average_enso").unlink()

    print("")
    print(" ===== Processing =====")

    if args.switch is None:
        print(" ON/OFF Switch : OFF")
    else:
        print(" ON/OFF Switch : ON")
        print(f" Assign to # of CONFS            : {args.switch}")
        anmr_enso['ONOFF'].fill(0)
        for i in range(len(args.switch)):
            anmr_enso['ONOFF'][args.switch[i]-1] = 1
    print(f" ON/OFF CONFS                    : {anmr_enso['ONOFF']}")

    if args.temp:
        print(" ON/OFF TEMP                     : ON ")
        print(f" TEMP is                         : {args.temp}")
        TEMP = args.temp  # temperature K
    else:
        print(" ON/OFF TEMP                     : OFF ")
        print(" TEMP is                         : 298.15 K [default]")
        TEMP = 298.15

    # For calculation the percentage of every CONFS
    anmr_enso = Boltzmann_enso(anmr_enso, TEMP)
    backup_enso = Boltzmann_enso(backup_enso, TEMP)

    names_anmr: list = list()
    if args.weights:
        # dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'), ('Energy', '<f8'),
        # ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8'), ('Total', '<f8'), ('Gibbs', '<f8'),
        # ('Qi', '<f8'), ('NEW_BW', '<f8')])

        # recalculation and copy to column 3 (BW)
        anmr_enso['BW'] = anmr_enso['NEW_BW']
        result_enso = np.copy(anmr_enso)

        # for normal condition, the ONOFF is 1(on) at all
        if result_enso.dtype.names:
            names_anmr = list(result_enso.dtype.names)
        result_enso['ONOFF'].fill(1)
        np.savetxt(args.file, result_enso[names_anmr[:8]], comments="", header="ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi",
                   fmt='%-6d %-4d %-4d %6.4f %11.7f %10.7f %10.7f %2.3f')

        print("")
        print("args.weights is ON and will only executive the percentage of each CONFS")
        print("Calculated and Weights of each CONFS from (Energy + mRRHO) ")
        print(f"And Save to the original file : {inFile}")
        print("Finished ...")
        print("")
        exit(0)

    # Average of CONFS
    # For after degeneracy and reduce the Gibbs free energy of ensemble
    avg_nums = int(np.sum(anmr_enso['ONOFF']))
    avg_fraction = 1/avg_nums

    print(f" the name of input file          : {inFile}")
    print(f" the name of input file energy   : {backupFile}")
    print(f" number of CONFS                 : {avg_nums:d}")
    print("")
    print(" ----- Average CONFS -----")

    if args.verbose:
        print("")
        print(" (1) Gibbs Free Energy of ensemble of average CONFS from Entropy ")
        print("")
        print(f"     Weight of every CONFS           : {avg_fraction: .4f}")
    avg_wlnw = math.log(avg_fraction)*(avg_fraction)

    if args.verbose:
        print(f"     Weight*ln(weight) for one CONF  : {avg_wlnw: .4f}")
        print(f"     Entropy of one CONF (cal/mol.K) : {avg_wlnw*(-Rcal): .4f}")  # nopep8
        print(f"     Temperature (K)                 : {TEMP: .2f}")
    Reduced_energy = avg_wlnw*(-Rcal)*TEMP/(-1000)

    if args.verbose:
        print(f"     Gibbs Free Energy  of one CONF  (kcal/mol)    : {Reduced_energy: .4f}")  # nopep8
        print(f"     Gibbs Free Energy  of all CONFS (kcal/mol)    : {Reduced_energy*avg_nums: .4f}")  # nopep8
        print("")
        print("")

    # For Before degeneracy and lift the Gibbs Free Energy of Ensemble
        print(" (2) Gibbs Free Energy of ensemble of average CONFS from Average Energy ")
    if args.verbose:
        print("")
        print(f"     ON/OFF of CONFS                 : {anmr_enso['ONOFF']}")
        print(f"     Eref of every CONFS (kcal/mol)  : {anmr_enso['Gibbs']}")

    insert_zero = False
    # enso_min_list = np.copy(anmr_enso['Gibbs']*anmr_enso['ONOFF'])

    for i in range(anmr_enso['Gibbs'].size):
        # is ON and Column 9 is delta Gibbs Free Energy
        if (anmr_enso['ONOFF'][i]) == 1 and (anmr_enso['Gibbs'][i]) == 0:
            insert_zero = True

    avg_Gibbs: npt.NDArray[np.float64] = anmr_enso['Gibbs'] * \
        anmr_enso['ONOFF']*(1/avg_nums)

    if args.verbose:
        print(f"     Weight of every CONFS           : {1/avg_nums}")
        print(f"     Eref*weight of every CONFS      : {avg_Gibbs}")
    Lift_energy: float = np.sum(avg_Gibbs).astype(float)

    if args.verbose:
        print(f"     Gibbs Free Energy (kcal/mol)    : {Lift_energy: .4f}\n")  # nopep8

    # Boltzmann of CONFS
    if args.verbose:
        print(" (3) Gibbs Free Energy of ensemble of Boltzmann of CONFS ")
    boltzmann_enso: npt.NDArray[np.float64] = (backup_enso['NEW_BW'] *
                                               anmr_enso['ONOFF']).astype(np.float64)
    sum_weight: float = np.sum(boltzmann_enso).astype(float)
    boltzmann_enso = boltzmann_enso/sum_weight

    # For Before degeneracy and lift the Gibbs free energy of ensemble
    if args.verbose:
        print("\n     Lift Gibbs Free Energy of ensemble using Boltzmann distribution ")
        print(f"     ON/OFF of CONFS                 : {anmr_enso['ONOFF']}")
        print(f"     Eref of every CONFS (kcal/mol)  : {backup_enso['Gibbs']}")
        print(f"     Weighting of each CONFS         : {anmr_enso['NEW_BW']}")

    Boltzmann_lift_energy: npt.NDArray[np.float64] = boltzmann_enso * \
        backup_enso['Gibbs']

    if args.verbose:
        print(
            f"     Eref*weight of every CONFS      : {Boltzmann_lift_energy}")
    Lift_boltzmann_energy: float = np.sum(Boltzmann_lift_energy).astype(float)
    if args.verbose:
        print(
            f"     Gibbs Free Energy  kcal/mol)    : {Lift_boltzmann_energy: .4f}")
        print("")

    print(f" (1) Gibbs Free Energy from Entropy (kcal/mol)                   : {Reduced_energy*avg_nums: .4f}")  # nopep8
    print(f" (2) Gibbs Free Energy from Average energy (kcal/mol)            : {Lift_energy: .4f}")  # nopep8
    print(f" (3) Gibbs Free Energy using Boltzmann distribution (kcal/mol)   : {Lift_boltzmann_energy: .4f}")  # nopep8
    print("")

    Gibbs_Free_Energy: float = Lift_energy - \
        Lift_boltzmann_energy + Reduced_energy*avg_nums
    if (Gibbs_Free_Energy) >= 0:
        rule = "(Forbidden)"
    else:
        rule = "(Allowed)"

    if args.verbose:
        print(" Total Gibbs Free Energy = G(Average energy) - G(Boltzmann distribution) + G(From Entropy)")
        print("                         = (2) - (3) + (1)                                                ")
    print(f" Total Gibbs Free Energy of ensemble of all CONFS (kcal/mol) : {Gibbs_Free_Energy: .4f}        ")  # nopep8
    print(f" Total Gibbs Free Energy of ensemble of all CONFS (Eh)       : {Gibbs_Free_Energy/Eh: .8f}    {rule}")  # nopep8

    if (Gibbs_Free_Energy <= 0):

        Reduced_energy_Eh = Reduced_energy / Eh
        # Gibbs_min is lowest energy of Gibbs Free Energy
        Gibbs_min = anmr_enso['Total'].min()
        Gibbs_Eh: npt.NDArray[np.float64] = np.copy(
            (anmr_enso['Total']-Gibbs_min)*anmr_enso['ONOFF'])
        insert_zero = False
        for i in range(anmr_enso['Total'].size):
            if (anmr_enso['ONOFF'][i]) == 1 and (Gibbs_Eh[i]) == 0.0:
                insert_zero = True
        Gibbs_Eh = Gibbs_Eh[Gibbs_Eh != 0]
        if insert_zero:
            Gibbs_Eh = np.insert(Gibbs_Eh, 0, 0)
        avg_Gibbs_Eh: float = np.average(Gibbs_Eh).astype(float)

        result_enso: npt.NDArray = np.copy(anmr_enso)

        # Conservation of Gibbs Free Energy, before ensemble after ensemble is the same Energy (not including the Reduced_Energy from Entropy)
        # dtype=[('ONOFF', '<i8'), ('NMR', '<i8'), ('CONF', '<i8'), ('BW', '<f8'), ('Energy', '<f8'),
        # ('Gsolv', '<f8'), ('mRRHO', '<f8'), ('gi', '<f8'), ('Total', '<f8'), ('Gibbs', '<f8'),
        # ('Qi', '<f8'), ('NEW_BW', '<f8')])
        result_enso['Energy'] = result_enso['Energy'] + \
            (avg_Gibbs_Eh-result_enso['Gibbs']/Eh)*(result_enso['ONOFF'])
        result_enso['mRRHO'] = result_enso['mRRHO'] + \
            Reduced_energy_Eh*result_enso['ONOFF']
        result_enso['BW'] = avg_fraction
        if result_enso.dtype.names:
            names_anmr = list(result_enso.dtype.names)
        np.savetxt(args.out, result_enso[names_anmr[:8]], comments="", header="ONOFF NMR  CONF BW      Energy        Gsolv      mRRHO      gi",
                   fmt='%-6d %-4d %-4d %6.4f %11.7f %10.7f %10.7f %2.3f')

        if args.verbose:
            print(
                "\n Average_enso is represented as like ensemble (Favor for reduced entropy of ensemble)")
            print(
                " the electronic energy and Entropy of every CONF is represented as single structure")


if __name__ == "__main__":
    main()

# test
# python3 ensoAnalysis.py -i tests/data/anmr_enso --new
# python3 ensoAnalysis.py -i tests/data/anmr_enso
#
