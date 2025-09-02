#!/usr/bin/env python
# import filecmp
import pytest
from censo_ext.Tools.utility import delete_all_files
from pathlib import Path
import argparse
import censo_ext.anmr as anmr
Dir = "tests/data/34.Ergocalciferol/04.Hydrogen"
outFile = "output.dat"
compare = "tests/compare/anmr_peaks.json"


def test_anmr_miss_args() -> None:
    x: dict = {}

    with pytest.raises(SystemExit) as e:
        anmr.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


x: dict = {"auto": True, "dir": Dir, "bobyqa": False, "mf": 500,
           "lw": None, "thr": None, "thrab": 0.025, "tb": 4, "mss": 9,
           "cutoff": 0.001, "show": False, "start": None, "end": None, "out": outFile}


def test_anmr_from_raw_data() -> None:
    x['average'] = False
    x['json'] = None
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    # assert filecmp.cmp(Dir/Path("peaks.json"), compare)


def test_anmr_average_on_json_off() -> None:
    x['average'] = True
    x['json'] = None
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    # delete_all_files("tests/data/34.Ergocalciferol/04.Hydrogen/peaks.json")   # Normal is necessary to remove the peaks.json but next method need this file


def test_anmr_average_on_json_on() -> None:
    x['average'] = True
    x['json'] = [-1]
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    delete_all_files(Dir / Path("peaks.json"), Dir/Path(outFile))


# if __name__ == "__main__":
#    import cProfile
#    cProfile.run("test_anmr_from_raw_data()", sort="cumtime")
