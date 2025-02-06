#!/usr/bin/env python3
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.molManipulate as molManipulate


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
        assert filecmp.cmp(dirName, compare) == True

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
    assert filecmp.cmp(args.out, compare) == True

    import shutil
    shutil.rmtree(Dir)


def test_molManipulate_multi():
    x: dict = {"separate": "tests/data/crest_conformers.xyz",
               "merge": None, "out": "tests/compare/output.xyz"}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        molManipulate.main(args)
    assert e.type == SystemExit
    assert e.value.code == 0   # for argprarse wrong


def test_molManipulate_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        molManipulate.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error
