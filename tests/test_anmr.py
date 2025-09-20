#!/usr/bin/env python
# import filecmp
import pytest
from censo_ext.Tools.utility import delete_all_files
from pathlib import Path
import argparse
import censo_ext.anmr as anmr
import filecmp
import platform
Dir_Ergo_H: Path = Path("tests/data/34.Ergocalciferol/04.Hydrogen")
Dir_Ergo_C: Path = Path("tests/data/34.Ergocalciferol/07.Carbon")
Dir_EA_H: Path = Path("tests/data/06.EthylAcetate/03.Censo")
outFile: Path = Path("output.dat")
compare_Ergo_H: Path = Path("tests/compare/anmr_peaks_H.json")
compare_Ergo_C: Path = Path("tests/compare/anmr_peaks_C.json")
_system: str = platform.system()


def test_anmr_miss_args() -> None:
    x: dict = {}

    with pytest.raises(SystemExit) as e:
        anmr.main(argparse.Namespace(**x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


x: dict = {"auto": True, "bobyqa": False, "mf": 500, "verbose": False,
           "lw": None, "thr": None, "thrab": 0.025, "tb": 4, "mss": 9,
           "cutoff": 0.001, "show": False, "start": None, "end": None, "out": outFile}


def test_anmr_H_from_raw_data() -> None:
    x['average'] = False
    x['json'] = None
    x['dir'] = Dir_Ergo_H
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    if _system == "Linux":
        compare: Path = Path("tests/compare/anmr_peaks_Ergo_H_ubuntu.json")
    elif _system == "Darwin":
        compare: Path = Path("tests/compare/anmr_peaks_Ergo_H_Darwin.json")

    import json
    import numpy as np
    with open(Dir_Ergo_H / Path("peaks.json"), "r") as jsonFile:
        source = np.array(json.load(jsonFile)[0])
    with open(compare, "r") as jsonFile:  # type: ignore
        denstination = np.array(json.load(jsonFile)[0])
    np.testing.assert_allclose(
        source, denstination, rtol=1e-100, atol=0)


def test_anmr_H_average_on_json_off() -> None:
    x['average'] = True
    x['json'] = None
    x['dir'] = Dir_Ergo_H
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    # delete_all_files("tests/data/34.Ergocalciferol/04.Hydrogen/peaks.json")   # Normal is necessary to remove the peaks.json but next method need this file


def test_anmr_H_average_on_json_on() -> None:
    x['average'] = True
    x['json'] = [-1]
    x['dir'] = Dir_Ergo_H
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 77868)
    delete_all_files(Dir_Ergo_H / Path("peaks.json"),
                     Dir_Ergo_H/Path(outFile))


def test_anmr_H_from_raw_data_EA() -> None:
    # this is Ethyl Acetate
    x['average'] = False
    x['json'] = None
    x['dir'] = Dir_EA_H
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 48310)
    if _system == "Linux":  # Need 2 min
        compare: Path = Path("tests/compare/anmr_peaks_EA_H_ubuntu.json")

    elif _system == "Darwin":  # Need 5 min
        compare: Path = Path("tests/compare/anmr_peaks_EA_H_Darwin.json")
    assert filecmp.cmp(Dir_EA_H/Path("peaks.json"), compare)  # type: ignore
    delete_all_files(Dir_EA_H / Path("peaks.json"), Dir_EA_H/Path(outFile))


def test_anmr_H_average_on_json_off_EA() -> None:
    # this is Ethyl Acetate
    x['average'] = True
    x['json'] = None
    x['dir'] = Dir_EA_H
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 48310)


def test_anmr_H_average_on_json_on_EA() -> None:
    # this is Ethyl Acetate
    x['average'] = True
    x['json'] = [-1]
    x['dir'] = Dir_EA_H
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 48310)
    delete_all_files(Dir_EA_H / Path("peaks.json"), Dir_EA_H/Path(outFile))


def test_anmr_C_from_raw_data() -> None:
    x['average'] = False
    x['json'] = None
    x['dir'] = Dir_Ergo_C
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 86411)
    if _system == "Linux":
        compare: Path = Path("tests/compare/anmr_peaks_Ergo_C_ubuntu.json")
    elif _system == "Darwin":
        compare: Path = Path("tests/compare/anmr_peaks_Ergo_C_Darwin.json")

    import json
    import numpy as np
    with open(Dir_Ergo_C / Path("peaks.json"), "r") as jsonFile:
        source = np.array(json.load(jsonFile)[0])
    with open(compare, "r") as jsonFile:  # type: ignore
        denstination = np.array(json.load(jsonFile)[0])
    np.testing.assert_allclose(
        source, denstination, rtol=1e-100, atol=0)


def test_anmr_C_average_on_json_off() -> None:
    x['average'] = True
    x['json'] = None
    x['dir'] = Dir_Ergo_C
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 86411)
    # delete_all_files("tests/data/34.Ergocalciferol/07.Carbon/peaks.json")   # Normal is necessary to remove the peaks.json but next method need this file


def test_anmr_C_average_on_json_on() -> None:
    x['average'] = True
    x['json'] = [-1]
    x['dir'] = Dir_Ergo_C
    assert anmr.main(argparse.Namespace(**x)).shape == (2, 86411)
    delete_all_files(Dir_Ergo_C / Path("peaks.json"), Dir_Ergo_C/Path(outFile))


if __name__ == "__main__":
    import cProfile
    cProfile.run("test_anmr_H_from_raw_data()", sort="cumtime")
