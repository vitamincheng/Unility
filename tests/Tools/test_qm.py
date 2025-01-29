#!/usr/bin/env python3
import pathlib
import argparse
from pathlib import Path
import pytest
from icecream import ic
from censo_ext.Tools.qm import *


def test_qm():

    x = {"out": "output.dat", "start": -
         0.5, "end": 10.5, "lw": 1, "mf": 500.0, "cutoff": 0.001, "debug": False}
    v = [964, 2775.76, 2768.20, 928]

    J = np.array([[0.0,   0.0,   0.0,   0.0],
                  [0.0,   0.0, 16.97,   0.0],
                  [0.0, 16.97,   0.0,   7.0],
                  [0.0,   0.0,   7.0,   0.0]])

    R_peak = qm_partial(v=v, J=J, nIntergals=1, idx0_nspins=1,
                        args=argparse.Namespace(**x))

    assert len(R_peak) == 16
    assert R_peak[0][0] == pytest.approx(np.float64(2751.6221950398317))

    assert R_peak[0][1] == pytest.approx(np.float64(0.004640626007360023))
    assert R_peak[len(R_peak) -
                  1][0] == pytest.approx(np.float64(2773.4946427349055))
    assert R_peak[len(R_peak) -
                  1][1] == pytest.approx(np.float64(0.09446798980250143))
