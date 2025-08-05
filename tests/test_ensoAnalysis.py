#!/usr/bin/env python3
import os
import pytest
import argparse
import censo_ext.ensoAnalysis as ensoAnalysis
import filecmp

file_anmr = "tests/data/34.Ergocalciferol/04.Hydrogen/anmr_enso"


def test_ensoAnalysis_miss_args():
    x: dict = {}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        ensoAnalysis.main(args)
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_ensoAnalysis_Hydrogen_miss():
    x: dict = {"file": file_anmr, "new": None}
    args = argparse.Namespace(**x)
    with pytest.raises(ValueError) as e:
        ensoAnalysis.main(args)
    assert str(e.value) == " Input file is exists but backup file is not exist. "


def test_ensoAnalysis_Hydrogen_miss_file():
    x: dict = {
        "file": "tests/data/34.Ergocalciferol/04.Hydrogen/anmr_enso1", "new": None}
    args = argparse.Namespace(**x)
    with pytest.raises(FileNotFoundError) as e:
        ensoAnalysis.main(args)
    assert str(
        e.value) == "tests/data/34.Ergocalciferol/04.Hydrogen/anmr_enso1 was not found or is a directory"


def test_ensoAnalysis_Hydrogen_new_read():
    x: dict = {"file": file_anmr, "new": True}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        ensoAnalysis.main(args)
    assert e.type is SystemExit
    assert e.value.code == 0  # for argparse error

    out_print = "result.txt"
    import sys
    with open(out_print, "w") as f:
        sys.stdout = f
        x: dict = {"file": file_anmr, "new": None, "switch": [*range(1, 55)],
                   "temp": 298.15, "weights": False, "complete": False}
        args = argparse.Namespace(**x)
        ensoAnalysis.main(args)
    sys.stdout = sys.__stdout__
    with open(out_print, "r") as f:
        lines = f.readlines()
    assert float(lines[-2].split()[-1]) == -0.1968
    assert lines[-1].split()[-1] == "(Allowed)"
    assert float(lines[-1].split()[-2]) == -0.00031361
    os.remove(out_print)
    os.remove(args.file+".backup")

    out_enso = "average_enso"
    compare = "tests/compare/test_average_enso"
    assert filecmp.cmp(out_enso, compare)
    os.remove(out_enso)


def test_ensoAnalysis_Hydrogen_new_read_miss_args():
    x: dict = {"file": file_anmr, "new": True}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        ensoAnalysis.main(args)
    assert e.type is SystemExit
    assert e.value.code == 0  # for argparse error

    out_print = "result.txt"
    import sys
    with open(out_print, "w") as f:
        sys.stdout = f
        x: dict = {"file": file_anmr, "new": None, "switch": None,
                   "temp": None, "weights": False, "complete": False}
        args = argparse.Namespace(**x)
        ensoAnalysis.main(args)
    sys.stdout = sys.__stdout__
    with open(out_print, "r") as f:
        lines = f.readlines()
    assert float(lines[-2].split()[-1]) == -0.1968
    assert lines[-1].split()[-1] == "(Allowed)"
    assert float(lines[-1].split()[-2]) == -0.00031361
    os.remove(out_print)
    os.remove(args.file+".backup")

    out_enso = "average_enso"
    compare = "tests/compare/test_average_enso"
    assert filecmp.cmp(out_enso, compare)
    os.remove(out_enso)


def test_ensoAnalysis_Hydrogen_new_read_complete():
    x: dict = {"file": file_anmr, "new": True}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        ensoAnalysis.main(args)
    assert e.type is SystemExit
    assert e.value.code == 0  # for argparse error

    out_print = "result.txt"
    import sys
    with open(out_print, "w") as f:
        sys.stdout = f
        x: dict = {"file": file_anmr, "new": None, "switch": None,
                   "temp": None, "weights": False, "complete": True}
        args = argparse.Namespace(**x)
        ensoAnalysis.main(args)
    sys.stdout = sys.__stdout__
    with open(out_print, "r") as f:
        lines = f.readlines()
    assert float(lines[-5].split()[-1]) == -0.1968
    assert lines[-4].split()[-1] == "(Allowed)"
    assert float(lines[-4].split()[-2]) == -0.00031361
    os.remove(out_print)
    os.remove(args.file+".backup")

    out_enso = "average_enso"
    compare = "tests/compare/test_average_enso"
    assert filecmp.cmp(out_enso, compare)
    os.remove(out_enso)


def test_ensoAnalysis_Hydrogen_new_read_weights():
    file_anmr_weights = "tests/data/34.Ergocalciferol/04.Hydrogen/anmr_enso_weights"
    x: dict = {"file": file_anmr_weights, "new": True}
    args = argparse.Namespace(**x)
    with pytest.raises(SystemExit) as e:
        ensoAnalysis.main(args)
    assert e.type is SystemExit
    assert e.value.code == 0  # for argparse error

    x: dict = {"file": file_anmr_weights, "new": None, "switch": None,
               "temp": None, "weights": True, "complete": True}
    args = argparse.Namespace(**x)

    with pytest.raises(SystemExit) as e:
        ensoAnalysis.main(args)
    assert e.type is SystemExit
    assert e.value.code == 0  # for argparse error
    os.remove(args.file+".backup")
