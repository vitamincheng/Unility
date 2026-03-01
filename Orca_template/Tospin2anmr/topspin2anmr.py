#!/usr/bin/env python3
from censo_ext.Tools.utility import save_simulation_spectra_file
import argparse
import numpy as np
from sys import argv as sysargv

descr = """
________________________________________________________________________________
| Input  : 1r file and procs (from Topspin)                                   
| Output : 1r.dat (Read by nmrplot.py)                                      
| Need   : procs                                                                              
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

    return parser.parse_args()


def search_string_in_file(fileName, string_to_search):
    """Search for the given string in file and return lines containing that string,
    along with line numbers"""
    list_of_results = []
    with open(fileName, 'r') as f:
        for line in f:
            if string_to_search in line:
                list_of_results.append(line.rstrip())
    return list_of_results


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    if args == argparse.Namespace():
        args = cml()

    print(descr)  # Program description
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

    from censo_ext.Tools.utility import IsExist
    IsExist("procs")

    print("Reading the procs file ")
    match_lines = search_string_in_file("procs", "ABSF1")
    SF1 = float(match_lines[0][match_lines[0].find(" ")+1:])
    print(match_lines[0])

    match_lines = search_string_in_file("procs", "ABSF2")
    SF2 = float(match_lines[0][match_lines[0].find(" ")+1:])
    print(match_lines[0])

    match_lines = search_string_in_file("procs", "$SI")
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

    with open(args.file, "rb") as f:
        byte: bytes = f.read(4)
        while byte:
            # Do stuff with byte.
            long: int = int.from_bytes(byte, byteorder='little', signed=True)
            outData.append([SF1 - i*step, long])
            i = i + 1
            byte = f.read(4)

    save_simulation_spectra_file(args.out, np.array(outData[::-1]))

    print(f"Coversion to anmr file {args.out}")
    print("Finished ...")


if __name__ == "__main__":
    main()


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
