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

    return float((median_075 - median_025 + median)*3)


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
    return (median - start_mean*20/19+median)*multi


def find_nearest(x_in: list[float] | npt.NDArray[np.float64], value) -> tuple[float, int]:
    """Finds the nearest value in a list to a given value.

    Args:
        x_in(list[float]): List of floats.
        value(float): Value to find the nearest to.

    Returns:
        tuple[float, int]: Tuple containing the nearest value and its index.
    """
    array: npt.NDArray[np.float64] = np.asarray(x_in)
    idx0: int = (np.abs(array - value)).argmin()
    return float(array[idx0]), idx0


def Boltzmann_Weighting(electron_Energy: npt.NDArray[np.float64], TEMP: float) -> npt.NDArray[np.float64]:

    # the unit of electron_Energy is kcal/mol
    # the unit of TEMP is K
    from censo_ext.Tools.Parameter import FACTOR

    # Gibbs_min is lowest energy of Gibbs Free Energy
    Gibbs_min: np.float64 = electron_Energy.min()
    Gibbs: npt.NDArray[np.float64] = np.array(electron_Energy-Gibbs_min)

    # Qi (each CONFS)
    Qi: npt.NDArray[np.float64] = np.array(np.exp(-Gibbs/(TEMP*FACTOR)))

    # Qall is sum of Qi
    Qall: np.float64 = np.sum(Qi)

    return Qi/Qall
