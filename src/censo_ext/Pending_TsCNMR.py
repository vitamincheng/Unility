#!/usr/bin/env python
import argparse
import numpy as np
from sys import argv as sysargv
from icecream import ic

descr = """
________________________________________________________________________________
|                                          [08.18.2024] vitamin.cheng@gmail.com
| Transform from ANMR .out file to .dat file (Modification)
| Usage    : TsCNMR.py [options]
| [options]
| Input    : -i input out file [default anmrc.out] 
| Output   : -o output dat file [default anmrc.dat]    
| Extra    : --extra xyz file for carbon's height 
|            (depend on different numbers of C-H and will show in nmrplot.py)  
| Debug    : -d debug mode and show the detail 
| Print    : -p print the final result on screen
| Packages : Tools
| Module   : anmrfile.py / ml4nmr.py / ASE library
|______________________________________________________________________________
"""


def cml(descr):
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--extra",
        dest="extra",
        action="store",
        required=False,
        help="Provide the name of input xyz file for carbon's height",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out",
        action="store",
        default="anmrc.dat",
        required=False,
        help="Provide the name of output file [default anmrc.dat]",
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        default="anmrc.out",
        required=False,
        help="Provide the name of input file [default anmrc.out]",
    )

    parser.add_argument(
        "-p",
        "--print",
        dest="print",
        action="store_true",
        help="Print the final result of calculation on screen",
    )

    parser.add_argument(
        "-d",
        "--debug",
        dest="debug",
        action="store_true",
        help="Debug mode or show the detail of calculation",
    )

    args = parser.parse_args()
    return args

########## GLOBAL DECLARATIONS ##########


Scan_width = 32 * 1024
Size_Frq = 238.8955
# LLimit_height = 0.1
Width_pp = Size_Frq / Scan_width
TMS_ppm = 14.5154
TMS_points = TMS_ppm / Width_pp
Line_width = 1
FWHM = Line_width / 1000
HMHM = FWHM / 2

########## END GLOBAL DECLARATIONS ##########


def Set_Peak_full_Spectra(x, multi, npData):
    for idx in range(Scan_width):
        npData[idx] = npData[idx] + \
            Lorentzian_Distribution(x-(idx*Width_pp-TMS_ppm), HMHM)*((multi+1))


def Lorentzian_Distribution(dx, mu):
    fx = 1 / 3.1415926 * (mu / (dx*dx + mu+mu))
    return fx


def main():
    args = cml(descr)
    if not args.print:
        print(descr)  # Program description
        print(f"    provided arguments: {" ".join(sysargv)}")

    from os.path import exists
    file_exists = exists(args.file)
    if not file_exists:
        print(f"    {args.file} , the file is not exist ...")
        print("  Exit and Close the program !!!")
        ic()
        exit(0)

    if not args.print:
        print("")
        print(f" Reading the {args.file} file ")

    import censo_ext.Tools.anmrfile as anmrfile
    anmrSJ = anmrfile.Anmr()
    anmrSJ.method_read_anmrSJ(args.file)

    if args.extra:
        import censo_ext.Tools.ml4nmr as ml4nmr
        mol, neighbors, bond_order = ml4nmr.read_mol_neighbors_bond_order(
            args.extra)
        List_nProton = [x for x in bond_order.values()]

    List_ppm = [ppm[3] for ppm in anmrSJ.anmrS]
    npData = np.zeros((Scan_width, 2))

    for idx in range(len(List_ppm)):
        if args.extra:
            Set_Peak_full_Spectra(
                List_ppm[idx], List_nProton[idx], npData.T[1])  # type: ignore
        else:
            Set_Peak_full_Spectra(List_ppm[idx], 1, npData.T[1])

    # npData.T[1] = npData.T[1]*0.1

    npData.T[0] = np.arange(Scan_width)*Width_pp-TMS_ppm

    threshold = 0.001
    outData = npData[np.logical_not(npData[:, 1] < threshold)]
    outData = np.insert(outData, 0, (npData[0][0], threshold), axis=0)
    outData = np.insert(outData, len(outData),
                        (npData[-1][0], threshold), axis=0)

    if args.out:
        np.savetxt(args.out, outData, fmt=' %2.6f %10.6e')
        print(f" Saved to the {args.out} file")

    if args.print:
        print(outData)


if __name__ == "__main__":
    main()
