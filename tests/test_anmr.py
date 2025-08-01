#!/usr/bin/env python3
import pytest
from censo_ext.Tools.utility import delete_all_files
from pathlib import Path
import argparse
import censo_ext.anmr as anmr
Dir = "tests/data/34.Ergocalciferol/04.Hydrogen"
outFile = "output.dat"


def test_anmr_miss_args() -> None:
    x: dict = {}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        anmr.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_anmr_from_raw_data() -> None:
    x: dict = {"auto": True, "average": False, "dir": Dir, "bobyqa": False, "mf": 500,
               "lw": None, "thr": None, "json": None, "thrab": 0.025, "tb": 4, "mss": 9,
               "cutoff": 0.001, "show": False, "start": None, "end": None, "out": outFile}
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)


def test_anmr_bobyqa_false_json_off() -> None:
    x: dict = {"auto": True, "average": True, "dir": Dir, "bobyqa": False, "mf": 500,
               "lw": None, "thr": None, "json": None, "thrab": 0.025, "tb": 4, "mss": 9,
               "cutoff": 0.001, "show": False, "start": None, "end": None, "out": outFile}
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    # delete_all_files("tests/data/34.Ergocalciferol/04.Hydrogen/peaks.json")   # Normal is necessary to remove the peaks.json but next method need this file


def test_anmr_bobyqa_true_json() -> None:
    x: dict = {"auto": True, "average": False, "dir": Dir, "bobyqa": True, "mf": 500,
               "lw": None, "thr": None, "json": [-1], "thrab": 0.025, "tb": 4, "mss": 9,
               "cutoff": 0.001, "show": False, "start": None, "end": None, "out": outFile}
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    delete_all_files(Dir / Path("peaks.json"), Dir/Path(outFile))


if __name__ == "__main__":
    import cProfile
    cProfile.run("test_anmr_from_raw_data()", sort="cumtime")
