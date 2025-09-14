#!/usr/bin/env python3
import sys
from os.path import exists
import argparse
import numpy as np
from sys import argv as sysargv

descr = """
________________________________________________________________________________
|                                          [01.06.2023] vitamin.cheng@gmail.com
| Input : 1r file (from Topspin)                                               
| Output : 1r.dat (Read by anmr program)                                       
|                                                                              
|______________________________________________________________________________
"""


def cml():
    # def cml():
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
        help="Provide name of the output file.[default 1r.dat]",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="1r",
        help="Provide input_file format [default 1r]",
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


print(descr)

args = cml()
print("    provided arguments: {}".format(" ".join(sysargv)))
print("")


file_exists = exists("procs")

if not file_exists:
    print("procs, the file is not exist ...")
    sys.exit()

print("Reading the procs file ")
match_lines = search_string_in_file("procs", "ABSF1")
SF1 = float(match_lines[0][match_lines[0].find(" ")+1:])
print(match_lines[0])

match_lines = search_string_in_file("procs", "ABSF2")
SF2 = float(match_lines[0][match_lines[0].find(" ")+1:])
print(match_lines[0])

# match_lines=search_string_in_file("procs", "$FTSIZE")
match_lines = search_string_in_file("procs", "$SI")
ftsize = float(match_lines[0][match_lines[0].find(" ")+1:])
print(match_lines[0])
print()

sw = SF1-SF2

step = sw/ftsize

print("Reading the 1r file ")
i = 0
outData = []
with open(args.file, "rb") as f:
    byte = f.read(4)
    while byte:
        # Do stuff with byte.
        long = int.from_bytes(byte, byteorder='little', signed=True)
        outData.append([SF1-i*step, long])
        i = i+1
        byte = f.read(4)


print("Coversion to anmr file (1r.dat)")

# data=np.array(outfile)

res = outData[::-1]

outfile_name = args.out
np.savetxt(outfile_name, res, fmt='%2.5f %12.5e')


# np.savetxt('1r.dat',data,fmt='%2.5e')


print("Finished ...")
