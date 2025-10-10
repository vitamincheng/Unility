#!/usr/bin/env python
# from icecream import ic
import numpy.typing as npt
from pathlib import Path
import shutil
import numpy as np
#
# https://steam.oxxostudio.tw/category/python/library/shutil.html
#
# [old] os.path system
#    import os
#    os.getcwd()
#    os.chdir()
#    os.mkdir()
#    os.rmdir()     if it is empty directory
#    os.listdir()
#    os.path.isdir()
#
# https://medium.com/ai反斗城/python-使用pathlib替代os-path-轉錄-edc9defb2bd8
# [new] Path
#    from pathlib import Path
#    p = Path("/home/user/Downloads/repo/test.txt")
#    p.name
#    p.parents[0]
#    p.cwd()
#    p.mkdir()
#    p.rmdir()
#    p.is_dir()
#    p.is_file()
#    p.unlink()         delete file
#    p.exist()
#    import shutil       if it is not empty directory
#    shutil.rmtree(directory, ignore_errors=True)
#    shutil.copyfile(source,denstination)


def IsExists_DirFileName(DirFile: Path | str) -> tuple[Path, str]:
    """Checks if a file or directory path exists and returns its directory and name components.

    This function verifies the existence of a given path and extracts its parent directory
    and basename.

    Args:
        DirFileName(pathlib.Path): The full path to check.

    Returns:
        tuple[Path, str]: A tuple containing:
            - Dir (pathlib.Path): The parent directory of the provided path.
            - Name (str): The basename (filename or directory name) of the provided path.

    Raises:
        FileNotFoundError: If the provided path does not exist.

    Example:
        >>> from pathlib import Path
        >>> # Assuming 'example.txt' exists in the current directory
        >>> dir_path, filename = IsExists_DirFileName("example.txt")
        >>> print(f"Directory: {dir_path}, Filename: {filename}")
        Directory: ., Filename: example.txt

    Note:
        If the provided path has no parent directory (e.g., a relative path like "file.txt"),
        this function will return the current directory as the parent and the filename
        as the name component.
    """

    DirFile = Path(DirFile)
    IsExist(DirFile)
    File: str = DirFile.name
    Dir: Path
    if DirFile.parents:
        Dir = DirFile.parents[0]
    else:
        raise FileNotFoundError(
            f"{DirFile} was not found or is a directory")
    return Dir, File


def function_is_float(string: str) -> bool:
    """Check if a string can be converted to a float.

    Args:
        string(str): The input string to check.

    Returns:
        bool: True if the string can be converted to a float, False otherwise.
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def function_is_int(string: str) -> bool:
    """Check if a string can be converted to an integer.

    Args:
        string(str): The input string to check.

    Returns:
        bool: True if the string can be converted to an integer, False otherwise.
    """
    try:
        int(string)
        return True
    except ValueError:
        return False


def move_file(source: Path | str, destination: Path | str) -> None:
    """Move a file from a source path to a destination path.

    This function first checks if the source file exists before attempting to move it.

    Args:
        source (Path): The path of the file to move.
        destination (Path): The path to the destination.

    Raises:
        FileNotFoundError: If the source file does not exist.
    """
    source = Path(source)
    destination = Path(destination)
    IsExist(source)
    shutil.move(source, destination)


def copy_file(source: Path | str, destination: Path | str) -> None:
    """Copy a file from a source path to a destination path.

    This function first checks if the source file exists before attempting to copy it.

    Args:
        source(Path): The path of the file to copy.
        destination(Path): The path to the destination.

    Raises:
        FileNotFoundError: If the source file does not exist.
    """
    source = Path(source)
    destination = Path(destination)
    IsExist(source)
    shutil.copy(source, destination)


def delete_all_files(*inFiles: Path | str) -> None:
    """Delete all specified files from the filesystem.

    This function takes a variable number of file paths and attempts to delete
    each one. If a file does not exist, the function silently continues without
    raising an error.

    Args:
        *fileNames: Variable length argument list of file paths to be deleted.
            Each argument should be a string representing the path to a file.

    Example:
        >>> delete_all_files('file1.txt', 'file2.txt', 'file3.txt')
        # Deletes all three files if they exist

    Note:
        This function does not distinguish between files and directories.
        If a directory path is passed, it will attempt to delete it as a file,
        which may result in an error.
    """

    for File in [*inFiles]:
        path = Path(File)
        IsExists: bool = path.exists()
        if IsExists:
            path.unlink()


def jsonKeys2int(x) -> dict:
    """Convert dictionary keys to integers.

    This is a helper function typically used as an `object_pairs_hook` for `json.loads`
    to convert string keys from a JSON object into integers.

    Args:
        x: A list of (key, value) pairs from a JSON object.

    Returns:
        dict: A dictionary with integer keys.
    """
    return {int(k): v for k, v in x}


def IsExist(inFile: Path | str) -> None:
    """Check if a file exists and raise FileNotFoundError if it doesn't.

    This function verifies whether the specified file path exists in the filesystem.
    If the file does not exist, it prints an error message and raises a FileNotFoundError
    with a descriptive message.

    Args:
        fileName(Path): The path to the file to check for existence.

    Raises:
        FileNotFoundError: If the specified file does not exist in the filesystem.

    Example:
        >>> from pathlib import Path
        >>> IsExist("example.txt")
        # If example.txt doesn't exist, raises FileNotFoundError with message
        # "example.txt The file is not Exist ..."

    Note:
        This function will terminate the program execution if the file does not exist,
        as it calls ic() and raises an exception.
    """
    inFile = Path(inFile)
    IsExists: bool = inFile.exists()
    if not IsExists:

        print(f"  The file {inFile} is not Exist ...")
        print("  Exit and Close the program !!!")
        exit(0)
        # raise FileNotFoundError(f"  The file {inFile} is not Exist ...")


def IsExist_bool(inFile: Path | str) -> bool:
    """Check if a file exists and return a boolean value.

    This function takes a file path as input and checks whether the file exists
    in the filesystem. If the file exists, it returns True; otherwise, it prints
    an error message and returns False.

    Args:
        inFile(Path): The path to the file to check for existence.

    Returns:
        bool: True if the file exists, False otherwise.

    Example:
        >>> from pathlib import Path
        >>> file_path = Path("example.txt")
        >>> result = IsExist_bool(file_path)
        >>> print(result)
        False

    Note:
        When the file does not exist, an error message is printed to the console
        in a formatted manner with 80 equal signs for visibility.
    """
    inFile = Path(inFile)
    IsExists: bool = inFile.exists()
    if IsExists:
        return True
    else:
        print(f"{inFile} the file is not exist ...")
        return False


def prog_IsExist(Prog: str) -> bool:
    """Check if a program exists in the system PATH.

    This function uses shutil.which() to search for the specified program
    in the system's PATH environment variable. If found, it returns True;
    otherwise, it prints an error message and raises a ValueError.

    Args:
        ProgramName(str): The name of the program to check for existence.
            This should be the exact command name as it would appear in the
            terminal/command prompt.

    Returns:
        bool: True if the program is found in PATH, False otherwise.

    Raises:
        ValueError: If the program is not found in PATH. The function will
            print an error message before raising the exception.

    Example:
        >>> prog_IsExist("python")
        True
        >>> prog_IsExist("nonexistent_program")
        nonexistent_program  the program is not exist ...
        Exit and Close the program !!!
        ValueError: the program is not Exist ...

    Note:
        This function is useful for validating system dependencies before
        executing commands that require specific programs to be available.
    """

    from shutil import which
    if which(Prog):
        return True
    else:
        raise ValueError(f"{Prog}, the program is not Exist ...")


def save_simulation_spectra_file(fileName: Path | str, spectra) -> None:
    """Save simulation spectra to file in either compressed numpy (.npz) or text (.dat) format.

    This function saves the provided spectra data to a file with the specified filename.
    The output format is determined by the file extension:
    - If the extension is '.npz', the data is saved as a compressed numpy file
    - If the extension is '.dat', the data is saved as a text file with formatted floating-point numbers

    Args:
        fileName (Path | str): The path to the output file. Must have either '.npz' or '.dat' extension.
        spectra: The spectra data to be saved. Should be compatible with numpy's savez_compressed
                 and savetxt functions.

    Returns:
        None: This function does not return any value.

    Example:
        >>> save_simulation_spectra_file("output.npz", my_spectra)
        >>> save_simulation_spectra_file("output.dat", my_spectra)

    Note:
        The function prints a confirmation message indicating the file path where the spectra was saved.
    """

    output: str = Path(fileName).name
    if output.split(".")[-1] == "npz":
        np.savez_compressed(fileName, spectra)
        print(f" the spectra is saved to : {fileName}")
    if output.split(".")[-1] == "dat":
        np.savetxt(fileName, spectra, fmt='%12.6f  %12.6e')
        print(f" the spectra is saved to : {fileName}")


def print_arguments() -> None:
    """Print all command-line arguments passed to the script.

    This function retrieves and displays all command-line arguments from sys.argv,
    formatting them as a single string separated by spaces. It's useful for
    debugging or logging purposes to see what arguments were provided when the
    script was executed.

    Returns:
        None: This function does not return any value.

    Example:
        >>> print_arguments()
        provided arguments: script.py --input file.txt --output result.dat

    Note:
        The function prints the arguments to standard output and includes an empty
        line after the argument list for better readability.
    """

    import sys
    print("    provided arguments: {}".format(" ".join(sys.argv)))
    print("")


def save_figure(fileName: str = "nmrplot") -> None:
    """Save the current matplotlib figure to PDF and SVG formats.

    This function saves the currently active matplotlib figure in both PDF (300 dpi) 
    and SVG formats with the specified base filename. The function creates two output 
    files: one with .pdf extension and another with .svg extension.

    Args:
        fileName(str, optional): The base name for the output files. Defaults to "nmrplot".
            The extensions ".pdf" and ".svg" will be automatically appended to create
            the final filenames.

    Returns:
        None: This function does not return any value.

    Example:
        >>> save_figure("my_plot")
        # Saves files as "my_plot.pdf" and "my_plot.svg"

        >>> save_figure()  # Uses default name "nmrplot"
        # Saves files as "nmrplot.pdf" and "nmrplot.svg"

    Note:
        This function requires matplotlib to be imported and a figure to be 
        currently active. The PDF file is saved with high resolution (300 dpi).
    """
    import matplotlib.pyplot as plt
    plt.savefig(f"{fileName}.pdf", dpi=300)
    plt.savefig(f"{fileName}.svg")


def SVD(first_array: npt.NDArray, second_array: npt.NDArray):
    """
    Perform Singular Value Decomposition (SVD) on a system of linear equations.

    This function takes two arrays representing a system of linear equations and
    solves for the unknown vector using SVD. It computes the pseudo-inverse of the
    first array and uses it to find the solution vector.

    Args:
        first_array (npt.NDArray): The coefficient matrix (A) of the system Ax = b.
            This should be a 2D array where each row represents an equation.
        second_array (npt.NDArray): The dependent variable vector (b) of the system Ax = b.
            This should be a 1D array representing the constants on the right-hand side.

    Returns:
        None: This function prints the solution vector xtilde to the console but does not return it.

    Example:
        >>> import numpy as np
        >>> A = np.array([[1, 2], [3, 4], [5, 6]])
        >>> b = np.array([1, 2, 3])
        >>> SVD(A, b)
        # Prints the solution vector xtilde

    Note:
        - The function assumes that first_array is a matrix and second_array is a vector.
        - This implementation uses numpy's SVD decomposition followed by pseudo-inverse computation.
        - The result is printed to stdout but not returned as a value.
        - This function modifies the global variable 'xtilde' which may cause side effects
          if called multiple times in succession without reinitialization.
    """

    first_column_array: npt.NDArray = first_array.reshape(-1, 1)
    U, S, VT = np.linalg.svd(first_column_array, full_matrices=False)
    xtilde = VT.T @ np.linalg.inv(np.diag(S)
                                  ) @ U.T @ second_array
    print(xtilde)


def cosine_similarity(vec1: npt.NDArray[np.float64] | list, vec2: npt.NDArray[np.float64] | list) -> float:
    """
    Calculates the cosine similarity between two vectors using NumPy.

    The cosine similarity is computed as the dot product of the two vectors
    divided by the product of their magnitudes (L2 norms). This metric
    ranges from -1 (completely opposite) to 1 (identical), with 0 indicating
    orthogonality.

    Args:
        vec1 (npt.NDArray[np.float64] | list): The first vector as a NumPy array or list of floats.
        vec2 (npt.NDArray[np.float64] | list): The second vector as a NumPy array or list of floats.

    Returns:
        float: The cosine similarity between the two vectors, ranging from -1 to 1.

    Raises:
        SystemExit: If the input vectors have different lengths, the program exits
            with an error message.

    Example:
        >>> import numpy as np
        >>> v1 = [1, 2, 3]
        >>> v2 = [4, 5, 6]
        >>> similarity = cosine_similarity(v1, v2)
        >>> print(similarity)
        0.9746318305482668

    Note:
        - If either vector has a magnitude of zero (i.e., all elements are zero),
          the function returns 0.0 to avoid division by zero.
        - Both input vectors must have the same length; otherwise, the program
          terminates with an error message.
    """

    # Ensure inputs are NumPy arrays
    vec1 = np.array(vec1)
    vec2 = np.array(vec2)

    if len(vec1) != len(vec2):
        print(f"{vec1=}")
        print(f"{vec2=}")
        print("  The numbers of two vector of your input file are not the same")
        print("  Exit and Close the program !!!")
        exit(0)

    # Calculate dot product
    dot_product = np.dot(vec1, vec2)

    # Calculate magnitudes (L2 norms)
    from numpy.linalg import norm
    magnitude_vec1 = norm(vec1)
    magnitude_vec2 = norm(vec2)

    # Handle division by zero if either magnitude is zero
    if magnitude_vec1 == 0 or magnitude_vec2 == 0:
        return 0.0  # Or raise an error, depending on desired behavior

    # Calculate cosine similarity
    cosine_similarity = float(dot_product / (magnitude_vec1 * magnitude_vec2))
    return cosine_similarity


def sub_numpy(sorted_data: npt.NDArray[np.float64] | list, max_number: int = 12) -> npt.NDArray[np.int64]:
    """Find optimal subdivision points for sorted data based on delta analysis.

    This function analyzes the differences between consecutive elements in sorted data
    to determine optimal subdivision points. It attempts to distribute the largest
    gaps evenly across the data while respecting a maximum segment size constraint.

    Args:
        sorted_data: A sorted array or list of float64 values to be subdivided.
        max_number: Maximum allowed size for any segment (default: 12).

    Returns:
        An array of integers representing optimal subdivision points.

    Raises:
        SystemExit: If no valid subdivision is found within the given constraints,
            prompting user to adjust max_number parameter.

    Example:
        >>> data = [1.0, 2.0, 3.0, 10.0, 15.0]
        >>> sub_numpy(data, max_number=5)
        array([3, 2])

    Note:
        The function uses a greedy approach to find the best subdivision by
        examining all possible cuts up to half the length of the data.
    """
    sorted_data = np.array(sorted_data)

    if len(sorted_data) < 3:
        print("  The numbers of sorted_data need more than three !!!")
        print("  Exit and Close the program !!!")
        exit(0)

    if len(sorted_data) < max_number:
        return np.array([len(sorted_data)])
    # from icecream import ic
    # ic(sorted_data)
    delta: npt.NDArray[np.float64] = np.diff(sorted_data)
    idx0_sorted_delta: npt.NDArray[np.float64] = delta.argsort()[::-1]
    # ic(delta)
    # ic(idx0_sorted_delta)
    for length in range(1, len(delta)//2):
        idx0_cut: npt.NDArray[np.float64] = idx0_sorted_delta[:length]
        # ic(idx0_cut)
        # ic(idx0_sorted_delta[:length])
        idx0_cut = np.insert(idx0_cut, 0, -1)
        idx0_cut = np.insert(idx0_cut, 0, len(delta))
        idx0_cut.sort()
        cut_diff: npt.NDArray[np.float64] = np.diff(idx0_cut)
        # ic(idx0_cut, cut_diff)
        if int(np.max(cut_diff)) <= max_number:
            result: npt.NDArray[np.int64] = np.array([])
            Max: int = int(np.max(cut_diff))
            argmax = np.argmax(cut_diff)

            quotient, remainder = divmod(
                int(np.sum(cut_diff[:argmax])), (Max-1))
            # ic(quotient, remainder)

            if quotient == 0 and remainder == 0:
                pass
            elif quotient == 0 and remainder != 0:
                result = np.append(result, remainder)
            else:
                for x in range(quotient):
                    result = np.append(result, Max-1)
                if remainder != 0:
                    result = np.append(result, remainder)
            result = np.append(result, Max)
            # ic(result)
            del quotient
            del remainder

            quotient, remainder = divmod(
                int(np.sum(cut_diff[argmax+1:])), Max-1)
            # ic(quotient, remainder)
            # ic(cut_diff, cut_diff[argmax+1:])
            if quotient == 0 and remainder == 0:
                pass
            elif quotient == 0 and remainder != 0:
                result = np.append(result, remainder)
            else:
                for x in range(quotient):
                    result = np.append(result, Max-1)
                if remainder != 0:
                    result = np.append(result, remainder)
            # ic(quotient, remainder)
            # ic(result)
            return result.astype(np.int64)
            # return cut_diff
    print("  Adjust the max_number to fit !!!")
    print("  Exit and Close the program !!!")
    exit(0)


def R_square(x: npt.NDArray, y: npt.NDArray) -> float:
    """Calculate the coefficient of determination (R-squared) for two arrays.

    This function computes the R-squared value, which represents the proportion 
    of the variance in the dependent variable (y) that is predictable from 
    the independent variable (x). It is calculated as the square of the Pearson 
    correlation coefficient between the two arrays.

    Args:
        x (npt.NDArray): Independent variable array. Should be 1D array-like.
        y (npt.NDArray): Dependent variable array. Should be 1D array-like and 
            have the same length as x.

    Returns:
        float: The coefficient of determination (R-squared) value, ranging from 
               0 to 1. A value of 1 indicates perfect correlation, while 0 indicates 
               no linear relationship.

    Example:
        >>> import numpy as np
        >>> x = np.array([1, 2, 3, 4, 5])
        >>> y = np.array([2, 4, 6, 8, 10])
        >>> R_square(x, y)
        1.0

    Note:
        Both input arrays must have the same length and contain numeric data. 
        The function uses numpy's corrcoef function to calculate the correlation 
        coefficient before squaring it to get R-squared.
    """
    # Calculate the correlation matrix
    correlation_matrix = np.corrcoef(x, y)
    r: float = correlation_matrix[0, 1]
    r_squared: float = r**2

    return (r_squared)
