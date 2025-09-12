#!/usr/bin/env python3
import sys
from os.path import exists
import argparse
import numpy as np
from sys import argv as sysargv
# from icecream import ic

descr = """
________________________________________________________________________________
|                                          [01.06.2023] vitamin.cheng@gmail.com
| Input  : 1r file and procs (from Topspin)                                   
| Output : 1r.dat                                       
| Need   : procs                                                                              
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
        default="1r",
        help="Provide input_file name (from topspin) [default 1r]",
    )

    args = parser.parse_args()
    return args


def search_string_in_file(fileName, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
    # nline = 0
    list_of_results = []
    with open(fileName, 'r') as read_obj:
        for line in read_obj:
            # nline += 1
            if string_to_search in line:
                list_of_results.append(line.rstrip())
    return list_of_results


print(descr)  # Program description
args = cml()
print("    provided arguments: {}".format(" ".join(sysargv)))
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

match_lines = search_string_in_file("procs", "$SI")
# match_lines=search_string_in_file("procs", "FTSIZE")
ftsize = float(match_lines[0][match_lines[0].find(" ")+1:])
print(match_lines[0])

# SF1    = 14.622
# SF2    = -2.6230
# ftsize = 512*1024


# print("Reading the procs file ")
# match_lines=search_string_in_file("procs", "SW_p")
# SW_p=float(match_lines[0][match_lines[0].find(" ")+1:])
# print(match_lines[0])
#
# match_lines=search_string_in_file("procs", "SF")
# SF=float(match_lines[0][match_lines[0].find(" ")+1:])
# print(match_lines[0])
#
# match_lines=search_string_in_file("procs", "$SI")
# ftsize=float(match_lines[0][match_lines[0].find(" ")+1:])
# print(match_lines[0])
#
#
# match_lines=search_string_in_file("procs", "$OFFSET")
# SF1=float(match_lines[0][match_lines[0].find(" ")+1:])
# print(match_lines[0])
#
# sw=SW_p/SF
#

sw = SF1-SF2

step = sw / ftsize
print(f"1 point of sw = {step}")
print()

print("Reading the 1r file ")
print()

i: int = 0
outData: list = []

infileName: str = args.file
outfileName: str = args.out

with open(infileName, "rb") as f:
    byte = f.read(4)
    while byte:
        # Do stuff with byte.
        long = int.from_bytes(byte, byteorder='little', signed=True)
        outData.append([SF1-i*step, long])
        i = i+1
        byte = f.read(4)

res: np.ndarray = np.array(outData[::-1])
np.savetxt(outfileName, res, fmt='%2.5f %12.5e')

print(f"Coversion to anmr file ({outfileName})")
print("Finished ...")

# use interpolation function returned by `interp1d`
# setting of plotting spectra
# Added start and end point

# start:float = -5.0
# end  :float = 15.0
# dpi  :int   = 10000
# res = np.insert(res,0,[start,0.0],axis=0)
# res = np.vstack((res,[end,0.0]))
#
# from scipy import interpolate
# f = interpolate.interp1d(res.T[0], res.T[1])
#
# xnew = np.linspace(start,end,int(end-start)*dpi+1)
# ynew = f(xnew)
#
# res_new = np.vstack((xnew,ynew))
#
# np.savetxt(outfile_name,res_new.T,fmt='%2.5f %12.5e')
#
# print("Coversion to anmr file (" + outfile_name + ")")
# print("Finished ...")
