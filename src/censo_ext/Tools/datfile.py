from __future__ import annotations
from scipy.interpolate import interp1d
from typing import Self
import numpy as np
import numpy.typing as npt
from pathlib import Path
from censo_ext.Tools.utility import IsExist
# from dataclasses import dataclass


class unit_conversion():

    def __init__(self, in_ppm: npt.NDArray[np.float64]) -> None:
        self.in_ppm: npt.NDArray[np.float64] = in_ppm
        self.args_ppm: dict = {idx: ppm for idx,
                               ppm in enumerate(in_ppm)}
        self.ppm_args: dict = {ppm: idx for idx,
                               ppm in enumerate(in_ppm)}

    def index(self, ppm) -> int:
        from censo_ext.Tools.spectra import find_nearest
        value, index = find_nearest(self.in_ppm, ppm)
        return index

    def ppm(self, index) -> float:
        return float(self.args_ppm[index])

    def ppm_scale(self) -> npt.NDArray[np.float64]:
        return self.in_ppm

    def ppm_limits(self) -> tuple[float, float]:
        return float(self.in_ppm.min()), float(self.in_ppm.max())


class Peaks_npz():
    """Class fo handling Peaks of NMR data files(npz format).
    """

    def __init__(self, uc: unit_conversion, file: Path | str = Path("peaks.npz")) -> None:
        # IsExist(file)
        self.__fileName: Path = Path(file)
        self.__peaks: npt.NDArray = np.array(
            [], dtype=[('cID', 'i8'), ('Start', 'f8'), ('End', 'f8'), ('Area', 'f8')])
        self.__uc: unit_conversion = uc

    def __len__(self):
        return len(self.__peaks)

    def method_ppm2cID(self, in_x_ppm: float):
        for peak in self.__peaks:
            if peak['Start'] > in_x_ppm > peak['End']:
                return peak['cID']
        return None

    def method_delete_cID(self, cIDs: list[int]) -> None:
        for cID in cIDs:
            if cID in self.__peaks['cID']:
                self.__peaks = self.__peaks[self.__peaks['cID'] != cID]
            else:
                print("  Delete element is wrong cID")
                print("  Exit and Close the program !!!")
                exit(0)

    def method_merge_cID(self, cIDs: list[int]) -> None:

        min_cID: int = np.array(cIDs).min()
        start, end = -99999, 99999
        for cID in sorted(cIDs):
            if cID in self.__peaks['cID']:
                args_x: npt.NDArray[np.intp] = np.argwhere(
                    self.__peaks['cID'] == cID)
                if start < self.__peaks[args_x[0]]['Start']:
                    start: float = self.__peaks[args_x[0]
                                                ]['Start'].max().astype(float)
                if end > self.__peaks[args_x[0]]['End']:
                    end: float = self.__peaks[args_x[0]
                                              ]['End'].min().astype(float)
            else:
                print("  Merge element is wrong cID")
                print("  Exit and Close the program !!!")
                exit(0)

        # assign the new cID of data
        args_x = np.argwhere(self.__peaks['cID'] == min_cID)
        Total_intensit = 0
        for cID in cIDs:
            Total_intensit += np.sum(
                self.__peaks[self.__peaks['cID'] == cID]['Area'])

        # Remove unnecessary entry
        for cID in cIDs:
            if (cID in self.__peaks['cID']):
                self.__peaks = self.__peaks[self.__peaks['cID'] != cID]

        self.__peaks = np.insert(
            self.__peaks, args_x[0], (min_cID, start, end, Total_intensit))

    def method_cut_cID(self, cID, intensit) -> None:
        # use ng.peakpick.pick from y_heighest 0.90 to down to two different peaks
        if cID in self.__peaks['cID']:
            args_x: npt.NDArray[np.intp] = np.argwhere(
                self.__peaks['cID'] == cID)
            l_peaks: float = self.__peaks[args_x][0]['Start'][0].astype(
                float)
            r_peaks: float = self.__peaks[args_x][0]['End'][0].astype(
                float)
            min: int = self.__uc.index(l_peaks)
            max: int = self.__uc.index(r_peaks)
            if min > max:
                min, max = max, min
            import nmrglue as ng
            y_highest = intensit[min:max+1].max()
            ratio: float = 0.90
            cut_peaks: npt.NDArray = np.array([])
            while (1):
                cut_peaks = ng.peakpick.pick(
                    data=intensit[min:max+1], pthres=y_highest*ratio, algorithm="downward")
                if (len(cut_peaks) >= 2):
                    break
                else:
                    ratio -= 0.10
            sorted_cut_peaks = np.sort(cut_peaks, order='VOL')
            start = int(sorted_cut_peaks['X_AXIS'][-1] + min)
            end = int(sorted_cut_peaks['X_AXIS'][-2] + min)
            if end < start:
                start, end = end, start
            cut_argmin: np.intp = np.argmin(intensit[start:end + 1])
            cut_center: float = self.__uc.ppm(start + cut_argmin)
            # remove the old entry and add two additional entry
            self.__peaks = self.__peaks[self.__peaks['cID'] != cID]

            self.__peaks = np.insert(
                self.__peaks, args_x[0], (cID, cut_center, r_peaks, intensit[min:start+cut_argmin].sum()))
            self.__peaks = np.insert(
                self.__peaks, args_x[0], (self.__peaks['cID'].max() + 1, l_peaks, cut_center, intensit[start+cut_argmin:max].sum()))
        else:
            print("  Merge element is wrong cID")
            print("  Exit and Close the program !!!")
            exit(1)

    def method_integrate(self, intensit) -> list[tuple[int, npt.NDArray, npt.NDArray]]:
        out_Data: list = []
        for cID, start, end, _ in self.__peaks:  # type: ignore
            min: int = self.__uc.index(start)
            max: int = self.__uc.index(end)
            if min > max:
                min, max = max, min

            # extract the peak
            peak_int: npt.NDArray = intensit[min:max + 1]
            peak_scale: npt.NDArray = self.__uc.ppm_scale()[min:max + 1]
            out_Data.append((cID, peak_int, peak_scale))
        return out_Data

    def method_load_Data(self, in_Data: list | npt.NDArray[np.float64]) -> None:
        self.__peaks = np.array(
            in_Data, dtype=[('cID', 'i8'), ('Start', 'f8'), ('End', 'f8'), ('Area', 'f8')])

    def method_read_file(self):
        file = Path(self.__fileName)
        from censo_ext.Tools.utility import IsExists_DirFileName
        _, Name = IsExists_DirFileName(file)
        self.__fileName = Path(Name)
        file_split: list[str] = Name.split(".")
        file_ext: str = file_split[1]

        if file_ext == "npz":
            in_Data = np.load(file)['arr_0']
            self.__peaks = in_Data
        else:
            print("  File extension is .npz file !!!")
            print("  Exit and Close the program !!!")
            exit(0)

    def get_cIDs_center_peaks(self) -> npt.NDArray[np.float64]:
        first: npt.NDArray[np.float64] = self.__peaks['cID']
        second: npt.NDArray[np.float64] = (
            self.__peaks['Start']+self.__peaks['End'])/2
        return np.stack((first, second))

    def get_peaks_integral_number(self) -> zip[tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]]:
        ppm: npt.NDArray[np.float64] = (
            self.__peaks['Start']+self.__peaks['End'])/2
        min_basic: np.float64 = np.min(self.__peaks['Area'])
        integral_number: npt.NDArray[np.float64] = self.__peaks['Area']/min_basic
        return zip(ppm, integral_number)

    def method_print(self) -> None:
        print(self.__fileName)
        print(self.__peaks)

    def method_save(self) -> None:
        from censo_ext.Tools.utility import save_simulation_spectra_file
        save_simulation_spectra_file(self.__fileName, self.__peaks)


class CensoDat():
    """Class for handling NMR data files (dat format).

    This class manages NMR data in dat format including reading, manipulating,
    and saving spectral data with various processing capabilities.
    """

    def __init__(self, file: Path | str = Path("anmr.dat")) -> None:
        """
        Initialize the CensoDat object.

        Args:
            file (Path | str): Path to the dat file. Defaults to Path("anmr.dat").

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
        # self.__dat [[chemical_shift, amplitude],[...]]

        IsExist(file)
        self.__fileName: Path = Path(file)
        from censo_ext.Tools.utility import IsExists_DirFileName
        _, fileName = IsExists_DirFileName(file)
        file_ext: str = fileName.split(".")[-1]

        if file_ext == "dat":
            self.__dat: npt.NDArray[np.float64] = np.genfromtxt(file)
        elif file_ext == "npz":
            in_Data = np.load(file)['arr_0']
            a, *b = in_Data.shape
            if a == 1:
                in_Data = in_Data[0]
            self.__dat: npt.NDArray[np.float64] = in_Data

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
            censoDat: Self = copy.deepcopy(self)
            censoDat.__dat[:, 1] = np.subtract(
                self.__dat[:, 1], other.__dat[:, 1])
        else:
            print("  Two dat file is not the same scale")
            print("  Exit and Close the program !!!")
            exit(0)
        return censoDat

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

        ext: str = Path(self.__fileName).name.split(".")[-1]
        if ext == "npz":
            np.savez_compressed(self.__fileName, self.__dat)
            print(f" the spectra is saved to : {self.__fileName}")
        if ext == "dat":
            np.savetxt(self.__fileName, self.__dat, fmt='%12.6f  %12.6e')
            print(f" the spectra is saved to : {self.__fileName}")

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
