#!/usr/bin/env python3
import re
from icecream import ic
import sys
from sys import argv as sysargv
import argparse
import numpy.typing as npt
import numpy as np

descr = """
________________________________________________________________________________
|                                          [01.06.2023] vitamin.cheng@gmail.com
| Input  : 1r.jdx file (JCAMP-DX format file)                                  
| Output : 1r.dat (Read by anmr program)                                       
|______________________________________________________________________________
"""


# def cml():
def cml():
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        #        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )  # argparse.RawDescriptionHelpFormatter) #,

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


print(descr)  # Program description
args = cml()
print(f"    provided arguments: {" ".join(sysargv)}")
print("")

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

if not args.file:
    print(f"{args.file}, the file is not exist ...")
    sys.exit(0)


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

outfile = []

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
    print("  Close and exit the program !!!")
    exit(0)

idx0_lines_start += 1
idx0_lines_end -= 1
outData: list = []
TotalNums: int = 0

for line in lines[idx0_lines_start:idx0_lines_end]:
    line = line.replace("-", " -")
    ChemShift: float = float(line.split()[0])
    nNums = len(line.split())-1
    for nums in range(nNums):
        # print(start-TotalNums*step, line.split()[nums+1])
        outData.append([float((start-TotalNums*step) / freq),
                       float(line.split()[nums+1])])
        TotalNums = TotalNums+1


print(f"Coversion to anmr file {args.out}")

res: npt.NDArray[np.float64] = np.array(outData[::-1])
np.savetxt(f"{args.out}", res, fmt='%2.5f %12.5e')

print("Finished ...")
