#!/usr/bin/env python3
from censo_ext.Tools.utility import save_simulation_spectra_file
import re
from icecream import ic
from sys import argv as sysargv
import argparse
import numpy as np

descr = """
________________________________________________________________________________
|                                          [01.06.2023] vitamin.cheng@gmail.com
| Input  : 1r.jdx file (JCAMP-DX format file)                                  
| Output : 1r.dat (Read by anmr program)                                       
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
        "-o",
        "--output",
        dest="out",
        action="store",
        required=False,
        default="1r.dat",
        help="Provide name of the output file without file ending. [default 1r.dat]",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="1r.jdx",
        help="Provide input_file name [default 1r.jdx]",
    )

    args: argparse.Namespace = parser.parse_args()
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
    print(f"    provided arguments: {" ".join(sysargv)}")

    # start  = 16.00485
    # end    = -3.990000
    # ftsize = 16384
    # step   = sw/ftsize
    # Topspin
    # start   Low field limit of spectrum (OFFSET)
    # end
    # sw      Spectral width (SW=SF1-SF2)
    # ftsize  size of real spectrum (SI)
    ################################################################

    from censo_ext.Tools.utility import IsExist
    IsExist(args.file)

    # for jcamp-dx : start of line
    print("Reading the jdx format file ")
    match_lines = search_string_in_file(args.file, "FIRSTX")
    start = float(match_lines[0][match_lines[0].find("=")+1:])
    print(match_lines[0])

    match_lines = search_string_in_file(args.file, "LASTX")
    end = float(match_lines[0][match_lines[0].find("=")+1:])
    print(match_lines[0])

    # match_lines = search_string_in_file(args.file, "DELTAX")
    # DELTAX = float(match_lines[0][match_lines[0].find("=")+1:])
    # print(match_lines[0])
    #
    match_lines = search_string_in_file(args.file, "FREQUENCY")
    freq = float(match_lines[0][match_lines[0].find("=")+1:])
    print(match_lines[0])

    match_lines = search_string_in_file(args.file, "NPOINTS")
    ftsize = float(match_lines[0][match_lines[0].find("=")+1:])
    print(match_lines[0])

    sw = start-end
    step = sw/ftsize
    ic(start, end)
    ic(sw, ftsize)
    ic(freq, step)

    lines: list[str] = open(args.file, "r").readlines()
    idx0_lines_start: int = 0
    idx0_lines_end: int = 0

    for idx0, line in enumerate(lines):
        if re.search(r"XYDATA", line):
            idx0_lines_start = idx0
        if re.search(r"END", line):
            idx0_lines_end = idx0

    if idx0_lines_start == 0 or idx0_lines_end == 0:
        print("  Your jdx format file have something wrong !!!")
        print("  Exit and Close the program !!!")
        exit(0)

    idx0_lines_start += 1
    idx0_lines_end -= 1
    outData: list = []
    TotalNums: int = 0

    for line in lines[idx0_lines_start:idx0_lines_end]:
        line: str = line.replace("-", " -")
        nNums: int = len(line.split())-1
        for nums in range(nNums):
            outData.append([float((start-TotalNums*step) / freq),
                           float(line.split()[nums+1])])
            TotalNums += 1

    save_simulation_spectra_file(args.out, np.array(outData[::-1]))

    print(f"Coversion to anmr file {args.out}")
    print("Finished ...")


if __name__ == "__main__":
    main()
