import pytest
import numpy as np
import numpy.typing as npt
from censo_ext.Tools.spectra import numpy_threshold, find_nearest, numpy_threshold_mean_3


in_np: npt.NDArray[np.float64] = np.array(
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 22, 39, 78, 12, 4, 4, 9, -10, -12])


def test_numpy_threshold():

    assert numpy_threshold(in_np, 10.0) == 246.31578947368422
    assert numpy_threshold(in_np, 3.0) == 73.89473684210526


def test_numpy_threshold_mean_3():

    assert numpy_threshold_mean_3(in_np) == 36


def test_find_nearest():

    assert (find_nearest(list(in_np), 24)) == (22, 19)
