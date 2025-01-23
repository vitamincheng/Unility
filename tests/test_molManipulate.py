import string
import pytest
from pathlib import Path
import argparse
import molManipulate


def test_molManipulate():
    x: dict = {"separate": "tests/data/crest_conformers3.xyz",
               "merge": None, "out": "tests/compare/output.xyz"}
    args = argparse.Namespace(**x)
    molManipulate.main(args)

    import os
    Dir: Path = Path("Separation")
    dirNames: list[str] = next(os.walk(Dir))[1]
    import filecmp
    for dirName in dirNames:
        compare = str("tests/compare/") + str(dirName)
        dcmp = filecmp.cmp(dirName, compare)
        assert (dcmp == True)

    x: dict = {"separate": None, "merge": ["Separation/1.xyz", "Separation/2.xyz"],
               "out": "Separation/1&2.xyz"}
    args = argparse.Namespace(**x)
    molManipulate.main(args)
    x: dict = {"separate": None, "merge": ["Separation/1&2.xyz", "Separation/3.xyz"],
               "out": "Separation/1&2&3.xyz"}
    args = argparse.Namespace(**x)
    molManipulate.main(args)

    import filecmp
    compare = "tests/compare/molManipulate.xyz"
    dcmp = filecmp.cmp(args.out, compare)
    assert (dcmp == True)

    import shutil
    shutil.rmtree(Dir)
