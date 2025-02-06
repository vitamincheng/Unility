import pytest
import numpy as np
from censo_ext.Tools.spectra import numpy_threshold_10, numpy_threshold_3, find_nearest, numpy_threshold_mean_3


in_np: np.ndarray = np.array(
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 22, 39, 78, 12, 4, 4, 9, -10, -12])


def test_numpy_threshold_10():

    assert numpy_threshold_10(in_np) == 246.31578947368422


def test_numpy_threshold_3():

    assert numpy_threshold_3(in_np) == 73.89473684210526


def test_numpy_threshold_mean_3():

    assert numpy_threshold_mean_3(in_np) == 36


def test_find_nearest():

    assert (find_nearest(in_np, 24)) == (22, 19)
