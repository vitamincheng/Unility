#!/usr/bin/env python3
import string
import pytest
from pathlib import Path
import argparse
import censo_ext.anmr as anmr
import filecmp
import os


def test_anmr():

    import censo_ext.anmr as anmr
    import argparse
    x: dict = {"auto": True, "average": False, "dir": "tests/data/04.Hydrogen",
               "bobyqa": True, "mf": 500, "thr": None, "json": [-1], "thrab": 0.025,
               "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    peaks = anmr.main(argparse.Namespace(**x))
    assert peaks.shape == (2, 77868)

    x: dict = {"auto": True, "average": True, "dir": "tests/data/04.Hydrogen",
               "bobyqa": False, "mf": 500, "thr": None, "json": None, "thrab": 0.025,
               "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    peaks = anmr.main(argparse.Namespace(**x))
    assert peaks.shape == (2, 77868)
