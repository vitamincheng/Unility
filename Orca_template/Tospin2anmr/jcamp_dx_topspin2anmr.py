#!/usr/bin/env python3
import sys
from os.path import exists
import argparse
import numpy as np
# from sys import argv as sysargv

descr = """
________________________________________________________________________________
|                                          [01.06.2023] vitamin.cheng@gmail.com
| Input  : 1r file and procs (from Topspin)                                   
| Output : 1r.dat (Read by anmr program)                                       
|                                                                              
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
        default="1r",
        help="Provide name of the output file without file ending.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        help="Provide input_file name (from topspin) ",
    )

    parser.add_argument(
        "-start",
        "--startppm",
        dest="start",
        action="store",
        required=False,
        default=0,
        type=float,
        help="Start plotting from x-axis '<start>' ppm.",
    )
    parser.add_argument(
        "-End",
        "--endppm",
        dest="end",
        action="store",
        required=False,
        default=11,
        type=float,
        help="End plotting at '<end>' ppm. Value of end has to be larger than "
             "value of start.",
    )
    parser.add_argument(
        "-ftsize",
        "--size",
        dest="ftsize",
        action="store",
        required=False,
        default=64*1024,
        type=int,
        help="Numbers of point in scannig size ",
    )

    args = parser.parse_args()
    return args


def search_string_in_file(file_name, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
    line_number = 0
    list_of_results = []
    with open(file_name, 'r') as read_obj:
        for line in read_obj:
            line_number += 1
            if string_to_search in line:
                list_of_results.append((line.rstrip()))
    return list_of_results


# def search_string_in_file(file_name, string_to_search):
#    """Search for the given string in file and return lines containing that string,
#    along with line numbers"""
#    line_number = 0
#    list_of_results = []
#    with open(file_name, 'r') as read_obj:
#        for line in read_obj:
#            line_number += 1
#            if string_to_search in line:
#                list_of_results.append((line_number, line.rstrip()))
#    return list_of_results

print(descr)  # Program description
# args, args_defaults = cml()
# print("    provided arguments: {}".format(" ".join(sysargv)))


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

file_exists = exists("procs")

if not file_exists:
    print("procs, the file is not exist ...")
    sys.exit()

# print("Reading the procs file ")
# match_lines=search_string_in_file("procs", "ABSF1")
# SF1=float(match_lines[0][match_lines[0].find(" ")+1:])
# print(match_lines[0])
#
# match_lines=search_string_in_file("procs", "ABSF2")
# SF2=float(match_lines[0][match_lines[0].find(" ")+1:])
# print(match_lines[0])
#
# match_lines=search_string_in_file("procs", "$SI")
# match_lines=search_string_in_file("procs", "FTSIZE")
# ftsize=float(match_lines[0][match_lines[0].find(" ")+1:])
# print(match_lines[0])


# sw   = SF1-SF2

# SF1    = 14.622
# SF2    = -2.6230
# ftsize = 512*1024


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

# for jcamp-dx :  end of line


step = sw / ftsize
print(f"1 point of sw = {step}")

print()


print("Reading the 1r file ")
print()

i = 0
outfile = []
infile_name = "1r"

outfile_name = infile_name+".dat"

with open(infile_name, "rb") as f:
    byte = f.read(4)
    while byte:
        # Do stuff with byte.
        long = int.from_bytes(byte, byteorder='little', signed=True)
        outfile.append([SF1-i*step, long])
        outfile.append([SF1-(i+0.33)*step, long])
        outfile.append([SF1-(i+0.66)*step, long])
        i = i+1
        byte = f.read(4)


print(f"Coversion to anmr file, {outfile_name}")

res = outfile[::-1]

np.savetxt(outfile_name, res, fmt='%2.5f %12.5e')
# np.savetxt(outfile_name,data,fmt='%2.5e')


print("Finished ...")
