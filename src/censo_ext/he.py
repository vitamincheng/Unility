#!/usr/bin/env python
from censo_ext.Tools.utility import print_arguments
import argparse
descr = """
________________________________________________________________________________
| For Generation of xyz molecule 
| Usages   : xyz.py <geometry> [options]
| [options]
|______________________________________________________________________________
"""


def cml() -> argparse.Namespace:
    """ Get args object from commandline interface. Needs argparse module."""
    parser = argparse.ArgumentParser(
        # description="descr",
        # # formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS)
    parser.add_argument(
        "-i",
        "--input",
        dest="file",
        action="store",
        required=False,
        default="traj.xyz",
        help="Provide one input xyz file [default traj.xyz]",
    )
    return parser.parse_args()


def main(args: argparse.Namespace = argparse.Namespace()) -> None:
    if args == argparse.Namespace():
        args = cml()
    print_arguments()


if __name__ == "__main__":
    main()
