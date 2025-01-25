#!/usr/bin/env python3
import numpy as np


def numpy_threshold_mean_3(x: np.ndarray) -> float:
    x = np.sort(x.flatten())
    median_025: float = x[int(len(x)*0.25)]
    median_075: float = x[int(len(x)*0.75)]
    median: float = x[int(len(x)*0.50)]
    thr: float = (median_075 - median_025 + median)*3
    return thr


def numpy_threshold_3(x: np.ndarray) -> float:
    x = np.sort(x.flatten())
    start_mean: float = float(np.mean(x[0:int(len(x)*0.05)]))
    median: float = x[int(len(x)*0.50)]
    thr: float = (median - start_mean*20/19+median)*3
    return thr


def numpy_threshold_10(x: np.ndarray) -> float:
    x = np.sort(x.flatten())
    start_mean: float = float(np.mean(x[0:int(len(x)*0.05)]))
    median: float = x[int(len(x)*0.50)]
    thr: float = (median - start_mean*20/19+median)*10
    return thr


def find_nearest(arrayIn, value):
    array: np.ndarray = np.asarray(arrayIn)
    idx0: int = (np.abs(array - value)).argmin()
    return array[idx0], idx0
