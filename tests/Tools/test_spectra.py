import pytest
import numpy as np
from censo_ext.Tools.spectra import numpy_threshold_10, numpy_threshold_3, find_nearest, numpy_threshold_mean_3


def test_numpy():

    a: np.ndarray = np.array(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 22, 39, 78, 12, 4, 4, 9, -10, -12])
    result = (numpy_threshold_10(a))
    assert result == 246.31578947368422
    result = (numpy_threshold_3(a))
    assert result == 73.89473684210526
    result = (numpy_threshold_mean_3(a))
    assert result == 36
    result = (find_nearest(a, 24))
    assert result == (22, 19)
