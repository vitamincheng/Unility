
from __future__ import annotations
from scipy.interpolate import interp1d
from typing import Self
import sys
import numpy as np
import numpy.typing as npt
from pathlib import Path

from censo_ext.Tools.utility import IsExist
# from dataclasses import dataclass


class CensoDat():
    """Class for handling NMR data files (dat format).

    This class manages NMR data in dat format including reading, manipulating,
    and saving spectral data with various processing capabilities.
    """

    def __init__(self, file: Path | str = Path("anmr.dat")) -> None:
        """
        Initialize the CensoDat object.

        Args:
            file (Path): Path to the dat file. Defaults to Path("anmr.dat").

        Raises:
            FileNotFoundError: If the specified file does not exist.
            ValueError: If the file is empty or contains invalid data.

        Example:
            >>> from censo_ext.Tools.datfile import CensoDat
            >>> dat = CensoDat()
            >>> dat = CensoDat("custom_file.dat")

        Note:
            The dat file should contain two columns: chemical_shift and amplitude.
            The data is loaded using numpy's genfromtxt function.
        """

        IsExist(file)
        self.__fileName: Path = Path(file)
        self.__dat: npt.NDArray[np.float64] = np.genfromtxt(file)
        # self.__dat [[chemical_shift, amplitude],[...]]

    def __len__(self) -> int:
        """
        Get the length of the data array.

        Returns:
            int: Number of data points.
        """
        return len(self.__dat)

    def __sub__(self, other: Self) -> CensoDat:
        """
        Subtract another CensoDat object from this one.

        Args:
            other (Self): Another CensoDat object to subtract.

        Returns:
            CensoDat: New CensoDat object with difference.
        """
        if np.array_equal(self.__dat[:, 0], other.__dat[:, 0]):
            import copy
            Res: Self = copy.deepcopy(self)
            Res.__dat[:, 1] = np.subtract(self.__dat[:, 1], other.__dat[:, 1])
        else:
            raise ValueError(" Two dat file is not the same scale")
        return Res

    def __repr__(self) -> str:
        """
        Return string representation of the data.

        Returns:
            str: Formatted string showing all data points.
        """
        Res: str = ""
        for x in self.__dat:
            Res = Res + f'{x[0]:>12.6f}  {x[1]:>12.6e}\n'
        return Res

    def method_save_dat(self) -> None:
        """Save the data to file.

        Writes the current data to the file specified in __fileName.

        Args:
            self: The instance of the class containing the data to be saved.

        Returns:
            None: This method does not return any value.

        Raises:
            IOError: If there is an issue opening or writing to the file.
            FileNotFoundError: If the specified file path does not exist.
        """

        with open(self.__fileName, "w") as f:
            sys.stdout = f
            print(self, end="")
        sys.stdout = sys.__stdout__

    def method_normalize_dat(self, start: float = -5.0, end: float = 15.0, dpi: int = 10000, highest: int = 10000) -> None:
        """Normalize the data to a specific range.

        This method normalizes the spectral data to a specified range by:
        1. Adjusting the range parameters based on the maximum ppm value
        2. Adding boundary points to the data
        3. Interpolating the data to the specified resolution
        4. Normalizing the y-values to the specified maximum

        Args:
            start: The start of the normalization range in ppm. Defaults to -5.0.
            end: The end of the normalization range in ppm. Defaults to 15.0.
            dpi: Data points per unit. Defaults to 10000.
            highest: Maximum value for normalization. Defaults to 10000.

        Example:
            >>> datfile.method_normalize_dat(start=-10, end=20, dpi=5000, highest=5000)
            Normalizes data from -10 to 20 ppm with 5000 dpi and maximum value of 5000

        Note:
            If the maximum ppm value exceeds 50 and end is less than 50,
            the method automatically adjusts parameters to -20 to 240 ppm with 500 dpi.
        """

        ppm_least: float = self.__dat[:, 0][-1]
        if ppm_least > 50 and end < 50:
            start, end, dpi = -20, 240, 500

        if len(self) != 0:

            res: npt.NDArray[np.float64] = self.__dat
            res = np.insert(res, 0, [start, 0.0], axis=0)
            res = np.vstack((res, [end, 0.0]))

            from scipy import interpolate
            f: interp1d = interpolate.interp1d(res[:, 0], res[:, 1])

            xnew: npt.NDArray[np.float64] = np.linspace(start, end, int(end-start)*dpi+1).astype(np.float64)  # nopep8
            ynew: npt.NDArray[np.float64] = f(xnew)

            res_new: npt.NDArray[np.float64] = np.vstack((xnew, ynew))
            res_new[1] = res_new[1] / np.max(res_new[1]) * highest
            self.__dat = res_new.T

    def set_fileName(self, file: Path | str) -> None:
        """
        Set the filename for output.

        Args:
            file (Path): New filename.
        """
        self.__fileName = Path(file)

    def get_fileName(self) -> Path:
        """
        Get the current filename.

        Returns:
            Path: Current filename.
        """
        return self.__fileName

    def get_Dat(self) -> npt.NDArray[np.float64]:
        """
        Get the raw data array.

        Returns:
            npt.NDArray: The data array.
        """
        return self.__dat
