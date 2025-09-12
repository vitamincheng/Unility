#!/usr/bin/env python
from icecream import ic
from pathlib import Path
import os
import shutil
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


def IsExists_DirFileName(DirFile: Path):
    """Checks if a file or directory path exists and returns its directory and name components.

    This function verifies the existence of a given path and extracts its parent directory
    and basename.

    Args:
        DirFileName (pathlib.Path): The full path to check.

    Returns:
        tuple[Path, str]: A tuple containing:
            - Dir (pathlib.Path): The parent directory of the provided path.
            - Name (str): The basename (filename or directory name) of the provided path.

    Raises:
        FileNotFoundError: If the provided path does not exist.

    Example:
        >>> from pathlib import Path
        >>> # Assuming 'example.txt' exists in the current directory
        >>> dir_path, filename = IsExists_DirFileName(Path("example.txt"))
        >>> print(f"Directory: {dir_path}, Filename: {filename}")
        Directory: ., Filename: example.txt

    Note:
        If the provided path has no parent directory (e.g., a relative path like "file.txt"),
        this function will return the current directory as the parent and the filename 
        as the name component.
    """

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
        string (str): The input string to check.

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
        string (str): The input string to check.

    Returns:
        bool: True if the string can be converted to an integer, False otherwise.
    """
    try:
        int(string)
        return True
    except ValueError:
        return False


def move_file(source: Path, destination: Path) -> None:
    """Move a file from a source path to a destination path.

    This function first checks if the source file exists before attempting to move it.

    Args:
        source (Path): The path of the file to move.
        destination (Path): The path to the destination.

    Raises:
        FileNotFoundError: If the source file does not exist.
    """
    IsExist(source)
    shutil.move(source, destination)


def copy_file(source: Path, destination: Path) -> None:
    """Copy a file from a source path to a destination path.

    This function first checks if the source file exists before attempting to copy it.

    Args:
        source (Path): The path of the file to copy.
        destination (Path): The path to the destination.

    Raises:
        FileNotFoundError: If the source file does not exist.
    """
    IsExist(source)
    shutil.copy(source, destination)


def delete_all_files(*inFiles) -> None:
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


def delete_file_bool(inFile: Path) -> bool:
    """Deletes a file if it exists and returns whether the operation was successful.

    This function checks if a file exists at the specified path, and if it does,
    removes it from the filesystem. It returns True if the file was successfully
    deleted, or False if the file did not exist or could not be deleted.

    Args:
        fileName (Path): The path to the file that should be deleted.

    Returns:
        bool: True if the file was successfully deleted, False if the file
            did not exist or could not be deleted.

    Example:
        >>> from pathlib import Path
        >>> file_path = Path("example.txt")
        >>> delete_file_bool(file_path)
        True
    """

    from os.path import exists
    IsExists: bool = exists(inFile)
    if IsExists:
        os.remove(inFile)
        return True
    else:
        return False


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


def save_dict_orcaS(inFile: Path, Data: dict) -> None:
    """Save a dictionary to a file in orcaS format.

    Each key-value pair is written on a new line, formatted as an integer key and a float value.

    Args:
        fileName (Path): The path to the output file.
        Data (dict): The dictionary to save.
    """
    with open(inFile, 'w') as f:
        for key, value in Data.items():
            f.write('%10d %12.5f \n' % (key, value))


def load_dict_orcaS(inFile: Path) -> dict:
    """Load a dictionary from a file in orcaS format.

    Each line is expected to contain an integer key and a float value, separated by whitespace.

    Args:
        fileName (Path): The path to the input file.

    Returns:
        dict: A dictionary with integer keys and float values.
    """
    IsExist(inFile)
    lines: list = open(inFile, "r").readlines()
    Data: dict[int, float] = {}
    for x in lines:
        Data[int(x.split()[0])] = float(x.split()[1])
    return Data


def IsExist(inFile: Path) -> None:
    """Check if a file exists and raise FileNotFoundError if it doesn't.

    This function verifies whether the specified file path exists in the filesystem.
    If the file does not exist, it prints an error message and raises a FileNotFoundError
    with a descriptive message.

    Args:
        fileName (Path): The path to the file to check for existence.

    Raises:
        FileNotFoundError: If the specified file does not exist in the filesystem.

    Example:
        >>> from pathlib import Path
        >>> IsExist(Path("example.txt"))
        # If example.txt doesn't exist, raises FileNotFoundError with message
        # "example.txt The file is not Exist ..."

    Note:
        This function will terminate the program execution if the file does not exist,
        as it calls ic() and raises an exception.
    """

    IsExists: bool = inFile.exists()
    if not IsExists:
        raise FileNotFoundError(f"  The file {inFile} is not Exist ...")


def IsExist_return_bool(inFile: Path) -> bool:
    """Check if a file exists and return a boolean value.

    This function takes a file path as input and checks whether the file exists
    in the filesystem. If the file exists, it returns True; otherwise, it prints
    an error message and returns False.

    Args:
        inFile (Path): The path to the file to check for existence.

    Returns:
        bool: True if the file exists, False otherwise.

    Example:
        >>> from pathlib import Path
        >>> file_path = Path("example.txt")
        >>> result = IsExist_return_bool(file_path)
        >>> print(result)
        False

    Note:
        When the file does not exist, an error message is printed to the console
        in a formatted manner with 80 equal signs for visibility.
    """
    IsExists: bool = inFile.exists()
    if IsExists:
        return True
    else:
        print("="*80)
        print(f"{inFile} the file is not exist ...")
        print("="*80)
        return False


def program_IsExist(Prog: str) -> bool:
    """Check if a program exists in the system PATH.

    This function uses shutil.which() to search for the specified program
    in the system's PATH environment variable. If found, it returns True;
    otherwise, it prints an error message and raises a ValueError.

    Args:
        ProgramName (str): The name of the program to check for existence.
            This should be the exact command name as it would appear in the
            terminal/command prompt.

    Returns:
        bool: True if the program is found in PATH, False otherwise.

    Raises:
        ValueError: If the program is not found in PATH. The function will
            print an error message before raising the exception.

    Example:
        >>> program_IsExist("python")
        True
        >>> program_IsExist("nonexistent_program")
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


def save_figure(fileName="nmrplot") -> None:
    """Save the current matplotlib figure to PDF and SVG formats.

    Args:
        fileName (str, optional): The base name for the output files.
            Defaults to "nmrplot". The extensions ".pdf" and ".svg" will be appended.
    """
    import matplotlib.pyplot as plt
    plt.savefig(f"{fileName}.pdf", dpi=300)
    plt.savefig(f"{fileName}.svg")
