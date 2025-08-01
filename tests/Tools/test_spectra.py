# import pytest
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.spectra import numpy_thr, find_nearest, numpy_thr_mean_3


inSample: npt.NDArray[np.float64] = np.array(
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 22, 39, 78, 12, 4, 4, 9, -10, -12])


def test_spectra_numpy_thr():

    assert numpy_thr(inSample, 10.0) == 246.31578947368422
    assert numpy_thr(inSample, 3.0) == 73.89473684210526


def test_spectra_numpy_thr_mean_3():

    assert numpy_thr_mean_3(inSample) == 36


def test_spectra_find_nearest():

    assert (find_nearest(list(inSample), 24)) == (22, 19)
