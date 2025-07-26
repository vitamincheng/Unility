#!/usr/bin/env python3
import numpy as np
import numpy.typing as npt


def numpy_threshold_mean_3(x_in: npt.NDArray[np.float64]) -> float:
    x: npt.NDArray[np.float64] = np.sort(x_in.flatten())
    median_025: float = x[int(len(x)*0.25)]
    median_075: float = x[int(len(x)*0.75)]
    median: float = x[int(len(x)*0.50)]
    thr: float = (median_075 - median_025 + median)*3
    return thr


def numpy_threshold(x_in: npt.NDArray[np.float64], multi: float) -> float:
    x: npt.NDArray[np.float64] = np.sort(x_in.flatten())
    start_mean: float = float(np.mean(x[0:int(len(x)*0.05)]))
    median: float = x[int(len(x)*0.50)]
    thr: float = (median - start_mean*20/19+median)*multi
    return thr


def find_nearest(List_In: list[float], value) -> tuple[float, int]:
    array: npt.NDArray[np.float64] = np.asarray(List_In)
    idx0: int = (np.abs(array - value)).argmin()
    return array[idx0], idx0
