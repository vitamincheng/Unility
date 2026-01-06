#!/usr/bin/env python
import pytest
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.qm import qm_partial, qm_full


def test_qm_miss_args() -> None:
    freq: list[float] = [964, 2775.76, 2768.20, 928, 120000]
    JCoups: npt.NDArray[np.float64] = np.array([[0.0,   0.0,   0.0,   0.0],
                                                [0.0,   0.0, 16.97,   0.0],
                                                [0.0, 16.97,   0.0,   7.0],
                                                [0.0,   0.0,   7.0,   0.0]])

    with pytest.raises(ValueError) as e:
        qm_partial(freq=freq, JCoups=JCoups, idx0_nspins=1,
                   _cutoff=0.001, _verbose=False)
    assert str(e.value) == "Your JCoups is Error"


def test_qm_partial_full() -> None:

    freq: list[float] = [964, 2775.76, 2768.20, 928]
    JCoups: npt.NDArray[np.float64] = np.array([[0.0,   0.0,   0.0,   0.0],
                                                [0.0,   0.0, 16.97,   0.0],
                                                [0.0, 16.97,   0.0,   7.0],
                                                [0.0,   0.0,   7.0,   0.0]])

    R_peak: list[tuple[float, float]] = qm_partial(
        freq=freq, JCoups=JCoups, idx0_nspins=1, _cutoff=0.001, _verbose=False)

    assert len(R_peak) == 16
    assert R_peak[0][0] == pytest.approx(2751.6221950398317)
    assert R_peak[0][1] == pytest.approx(0.03711639618338153)
    assert R_peak[-1][0] == pytest.approx(2773.4946427349055)
    assert R_peak[-1][1] == pytest.approx(0.7555686087601723)

    R_peak = qm_full(freq=freq, JCoups=JCoups, _verbose=False, _cutoff=0.001)

    assert len(R_peak) == 36
    assert R_peak[0][0] == pytest.approx(924.4933121601566)
    assert R_peak[0][1] == pytest.approx(0.9961959539774445)
    assert R_peak[-1][0] == pytest.approx(931.493373552777)
    assert R_peak[-1][1] == pytest.approx(1.003803901594175)
