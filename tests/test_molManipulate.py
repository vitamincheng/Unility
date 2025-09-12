#!/usr/bin/env python
import pytest
from pathlib import Path
import argparse
import censo_ext.molManipulate as molManipulate


def test_molManipulate():
    x: dict = {"separate": Path("tests/data/crest_conformers3.xyz"),
               "merge": None, "out": Path("tests/compare/output.xyz")}
    molManipulate.main(argparse.Namespace(**x))

    Dir: Path = Path("Separation")
    dirNames: list[str] = [x for x in Dir.walk()][0][1]
    import filecmp
    for dirName in dirNames:
        compare: Path = Path("tests/compare") / Path(dirName)
        assert filecmp.cmp(dirName, compare)

    x: dict = {"separate": None, "merge": ["Separation/1.xyz", "Separation/2.xyz"],
               "out": Path("Separation/1&2.xyz")}
    args = argparse.Namespace(**x)
    molManipulate.main(args)
    x: dict = {"separate": None, "merge": ["Separation/1&2.xyz", "Separation/3.xyz"],
               "out": Path("Separation/1&2&3.xyz")}
    args = argparse.Namespace(**x)
    molManipulate.main(args)

    import filecmp
    compare = Path("tests/compare/molManipulate.xyz")
    assert filecmp.cmp(args.out, compare)

    import shutil
    shutil.rmtree(Dir)


def test_molManipulate_multi():
    x: dict = {"separate": Path("tests/data/crest_conformers.xyz"),
               "merge": None, "out": Path("tests/compare/output.xyz")}

    with pytest.raises(SystemExit) as e:
        molManipulate.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 0   # for argprarse wrong


def test_molManipulate_miss_args():
    x: dict = {}
    with pytest.raises(SystemExit) as e:
        molManipulate.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error
