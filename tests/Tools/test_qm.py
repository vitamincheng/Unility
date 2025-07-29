#!/usr/bin/env python3
import argparse
from pathlib import Path
import pytest
import numpy as np
import numpy.typing as npt
from icecream import ic
from censo_ext.Tools.qm import qm_partial, qm_full


def test_qm_miss_args():
    x = {"out": "output.dat", "start": -
         0.5, "end": 10.5, "lw": 1, "mf": 500.0, "cutoff": 0.001, "debug": False}
    v: list[float] = [964, 2775.76, 2768.20, 928, 120000]

    J: npt.NDArray[np.float64] = np.array([[0.0,   0.0,   0.0,   0.0],
                                           [0.0,   0.0, 16.97,   0.0],
                                           [0.0, 16.97,   0.0,   7.0],
                                           [0.0,   0.0,   7.0,   0.0]])

    with pytest.raises(ValueError) as e:
        R_peak: list[tuple[float, float]] = qm_partial(v=v, J=J, nIntergals=1,
                                                       idx0_nspins=1, args=argparse.Namespace(**x))
    assert str(e.value) == "Your JCoup is Error"


def test_qm_partial_full():

    x = {"out": "output.dat", "start": -
         0.5, "end": 10.5, "lw": 1, "mf": 500.0, "cutoff": 0.001, "debug": False}
    v: list[float] = [964, 2775.76, 2768.20, 928]

    J: npt.NDArray[np.float64] = np.array([[0.0,   0.0,   0.0,   0.0],
                                           [0.0,   0.0, 16.97,   0.0],
                                           [0.0, 16.97,   0.0,   7.0],
                                           [0.0,   0.0,   7.0,   0.0]])

    R_peak: list[tuple[float, float]] = qm_partial(v=v, J=J, nIntergals=1, idx0_nspins=1,
                                                   args=argparse.Namespace(**x))

    assert len(R_peak) == 16
    assert R_peak[0][0] == pytest.approx(2751.6221950398317)
    assert R_peak[0][1] == pytest.approx(0.004640626007360023)
    assert R_peak[-1][0] == pytest.approx(2773.4946427349055)
    assert R_peak[-1][1] == pytest.approx(0.09446798980250143)

    x = {"out": "output.dat", "start": -
         0.5, "end": 10.5, "lw": 1, "mf": 500.0, "cutoff": 0.001, "debug": False, "bobyqa": True}
    R_peak = qm_full(v=v, J=J, nIntergals=1, args=argparse.Namespace(**x))

    assert len(R_peak) == 36
    assert R_peak[0][0] == pytest.approx(924.4933121601566)
    assert R_peak[0][1] == pytest.approx(0.03113112356179515)
    assert R_peak[-1][0] == pytest.approx(931.493373552777)
    assert R_peak[-1][1] == pytest.approx(0.03136887192481798)
