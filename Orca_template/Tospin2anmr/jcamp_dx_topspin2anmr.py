#!/usr/bin/env python3
from censo_ext.Tools.utility import save_simulation_spectra_file_dat_npz
from sys import argv as sysargv
import argparse
import numpy as np

descr = """
________________________________________________________________________________
|                                          [01.06.2023] vitamin.cheng@gmail.com
| Input  : 1r file and procs (from Topspin)                                   
| Output : 1r.dat (Read by anmr program)                                       
|                                                                              
|______________________________________________________________________________
"""


def cml():
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
        default="1r.dat",
        help="Provide name of the output file without file ending [default 1r.dat] ",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="1r",
        help="Provide input_file name (from topspin) [default 1r] ",
    )

    args = parser.parse_args()
    return args


def search_string_in_file(file_name, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
    list_of_results = []
    with open(file_name, 'r') as f:
        for line in f:
            if string_to_search in line:
                list_of_results.append((line.rstrip()))
    return list_of_results


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print(descr)  # Program description
    args = cml()
    print(f"    provided arguments: {" ".join(sysargv)}")
    print("")

    from censo_ext.Tools.utility import IsExist
    IsExist("procs")

    # for jcamp-dx : start of line
    print("Reading the procs file ")
    match_lines = search_string_in_file("procs", "SW_p")
    SW_p = float(match_lines[0][match_lines[0].find(" ")+1:])
    print(match_lines[0])

    match_lines = search_string_in_file("procs", "##$SF")
    SF = float(match_lines[0][match_lines[0].find(" ")+1:])
    print(match_lines[0])

    match_lines = search_string_in_file("procs", "$SI")
    ftsize = float(match_lines[0][match_lines[0].find(" ")+1:])
    print(match_lines[0])

    match_lines = search_string_in_file("procs", "$OFFSET")
    SF1 = float(match_lines[0][match_lines[0].find(" ")+1:])
    print(match_lines[0])

    sw = SW_p/SF

    step = sw / ftsize
    print(f"1 point of sw = {step}")
    print()

    print("Reading the 1r file ")
    print()

    i: int = 0
    outData: list = []

    with open(args.file, "rb") as f:
        byte: bytes = f.read(4)
        while byte:
            # Do stuff with byte.
            long: int = int.from_bytes(byte, byteorder='little', signed=True)
            outData.append([SF1 - i*step, long])
            outData.append([SF1 - (i+0.33)*step, long])
            outData.append([SF1 - (i+0.66)*step, long])
            i += 1
            byte = f.read(4)

    save_simulation_spectra_file_dat_npz(args.out, np.array(outData[::-1]))

    print(f"Coversion to anmr file, {args.out}")
    print("Finished ...")


if __name__ == "__main__":
    main()
