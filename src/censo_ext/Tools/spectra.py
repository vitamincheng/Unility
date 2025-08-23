#!/usr/bin/env python
import numpy as np
import numpy.typing as npt


def numpy_thr_mean_3(x_in: npt.NDArray[np.float64]) -> float:
    """Calculates a threshold based on the 25th, 27th and 28th percentiles.

    The threshold is computed as: (P75 - P25 + median) * 3 where P25, P75 are
    the 25th and 75th percentiles respectively.

    Args:
        x_in: Input array of floats to compute threshold from.

    Returns:
        Calculated threshold value.
    """

    x: npt.NDArray[np.float64] = np.sort(x_in.flatten())
    median_025: float = x[int(len(x)*0.25)]
    median_075: float = x[int(len(x)*0.75)]
    median: float = x[int(len(x)*0.50)]
    thr: float = (median_075 - median_025 + median)*3
    return thr


def numpy_thr(x_in: npt.NDArray[np.float64], multi: float) -> float:
    """Calculates a threshold based on the mean and median of the input array.

    The threshold is computed as: (median - start_mean*20/19 + median) * multiplier
    where start_mean is the mean of the first 5% of sorted data.

    Args:
        x_in: Input array of floats to compute threshold from.
        multi: Multiplier for the threshold calculation.

    Returns:
        Calculated threshold value.
    """

    x: npt.NDArray[np.float64] = np.sort(x_in.flatten())
    start_mean: float = float(np.mean(x[0:int(len(x)*0.05)]))
    median: float = x[int(len(x)*0.50)]
    thr: float = (median - start_mean*20/19+median)*multi
    return thr


def find_nearest(x_in: list[float], value) -> tuple[float, int]:
    """Finds the nearest value in a list to a given value.

    Args:
        x_in (list[float]): List of floats.
        value (float): Value to find the nearest to.

    Returns:
        tuple[float, int]: Tuple containing the nearest value and its index.
    """
    array: npt.NDArray[np.float64] = np.asarray(x_in)
    idx0: int = (np.abs(array - value)).argmin()
    return float(array[idx0]), idx0
