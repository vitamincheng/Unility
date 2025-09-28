#! /usr/bin/env python
import nmrglue as ng
from nmrglue.fileio.fileiobase import unit_conversion
import numpy as np
import numpy.typing as npt
import argparse
from censo_ext.Tools.utility import print_descr


descr = """
________________________________________________________________________________
| Usage: ng2npz.py <geometry> [options]
| Transform Bruker fid file to npz file using nmrglue module
| [Options]
| Input     : -i the pdata path(under 1r folder) [required]
| Output    : -o the output npz file [default output.npz]
| Start     : -start start point of chemical shift [default from data]
| End       : -end end point of chemical shift [default from data]
| phase     : -p --phase [default 1.0]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface.
        Needs argparse module."""
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
        add_help=False
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="path",
        action="store",
        type=str,
        required=False,
        help="Provide the path of your pdata (under 1r folder) ",
    )

    parser.add_argument(
        "-o",
        "--out",
        dest="out",
        action="store",
        type=str,
        required=False,
        default="output.npz",
        help="Provide the filename of output [default output.npz] ",
    )

    parser.add_argument(
        "-start",
        dest="start",
        action="store",
        type=float,
        required=False,
        default=None,
        help="start point of chemical shift [default from data]",
    )

    parser.add_argument(
        "-end",
        dest="end",
        action="store",
        type=float,
        required=False,
        default=None,
        help="end point of chemical shift [default from data]",
    )

    parser.add_argument(
        "-p",
        "--phase",
        dest="phase",
        action="store",
        type=float,
        required=False,
        default=1.0,
        help="phase of spectra (1 to -1) [default 1.0]",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def main(args: argparse.Namespace = argparse.Namespace()) -> None:

    # read in the Bruker data
    if args == argparse.Namespace():
        args = cml()
    print_descr(descr)

    if not args.path:
        args.path = "../../Simulation/bmse000510/nmr/set01/1H/pdata/1"
        # args.path = "../../Simulation/bmse000510/nmr/set01/13C/pdata/1"
        # args.path = "../../Simulation/bmse000510/nmr/set01/DEPT_90/pdata/1"
        # args.path = "../../Simulation/bmse000510/nmr/set01/DEPT_135/pdata/1"
    dic, data = ng.bruker.read_pdata(args.path)
    udic: dict = ng.bruker.guess_udic(dic, data)

    C = ng.convert.converter()
    C.from_bruker(dic, data, udic)
    pipe_fid_fileName: str = ".1d_pipe.fid"
    ng.pipe.write(pipe_fid_fileName, *C.to_pipe(), overwrite=True)
    dic, data = ng.pipe.read(pipe_fid_fileName)
    data = data.real*args.phase  # type: ignore
    uc: unit_conversion = ng.pipe.make_uc(dic, data)

    # end ---------+--------- start
    # args.end                args.start
    # limit_ppm[0]            limit_ppm[1]

    limits_ppm = uc.ppm_limits()
    ppm: npt.NDArray[np.float64] = np.linspace(
        limits_ppm[0], limits_ppm[1], data.shape[0])
    if args.start is None or args.end is None:
        args.end, args.start = uc.ppm_limits()

    output: npt.NDArray[np.float64] = np.vstack((ppm, np.real(data))).T[::-1]
    from censo_ext.Tools.utility import save_simulation_spectra_file
    save_simulation_spectra_file(args.out, output)
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(pipe_fid_fileName)


if __name__ == "__main__":
    main()
