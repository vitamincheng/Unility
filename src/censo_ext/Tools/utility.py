#!/usr/bin/env python3
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
#
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
#    p.exist()
#    import shutil       if it is not empty directory
#    shutil.rmtree(directory, ignore_errors=True)


def IsExists_DirFileName(DirFileName: Path):
    IsExist(DirFileName)
    Name: str = DirFileName.name
    Dir: Path
    if DirFileName.parents:
        Dir = DirFileName.parents[0]
    else:
        print(DirFileName, " the DirFileName is not exist ...")
        print(" Exit and close the program !!! ")
        ic()
        raise FileNotFoundError(
            str(DirFileName) + " was not found or is a directory")
    return Dir, Name


def function_is_float(string: str) -> bool:
    try:
        float(string)
        return True
    except ValueError:
        return False


def function_is_int(string: str) -> bool:
    try:
        int(string)
        return True
    except ValueError:
        return False


def move_file(source: Path, destination: Path) -> None:
    IsExist(source)
    shutil.move(source, destination)


def copy_file(source: Path, destination: Path) -> None:
    IsExist(source)
    shutil.copy(source, destination)


def delete_all_files(*fileNames) -> None:
    from os.path import exists
    for file in [*fileNames]:
        IsExists: bool = exists(file)
        if IsExists:
            os.remove(file)


def delete_file_bool(fileName: Path) -> bool:
    from os.path import exists
    IsExists: bool = exists(fileName)
    if IsExists:
        os.remove(fileName)
        return True
    else:
        return False


def jsonKeys2int(x) -> dict:
    return {int(k): v for k, v in x}


def save_dict(fileName: Path, Data: dict) -> None:
    with open(fileName, 'w') as f:
        for key, value in Data.items():
            f.write('%12.5f %12.5e\n' % (key, value))


# def save_BOBYQA_orcaS(fileName: Path, Data: dict) -> None:
#    with open(fileName, 'w') as f:
#
#        for key, value in Data.items():
#            if type(value) is list:
#                f.write('%10d %12.5f %10d \n' % (key, value[0], value[1]))
#            else:
#                f.write('%10d %12.5f %10d \n' % (key, value, 0))
#
#
# def load_BOBYQA_orcaS(fileName: Path) -> dict:
#    IsExist(fileName)
#    lines = open(fileName, "r").readlines()
#
#    Data: dict = {}
#    for x in lines:
#        if len(x.split()) == 2:
#            Data[int(x.split()[0])] = [float(x.split()[1])]
#        elif len(x.split()) == 3:
#            Data[int(x.split()[0])] = [float(x.split()[1]), int(x.split()[2])]
#        else:
#            print("something wrong in your orcaS.out")
#            print("Exit to the program")
#            ic()
#    return Data


def save_dict_orcaS(fileName: Path, Data: dict) -> None:
    with open(fileName, 'w') as f:
        for key, value in Data.items():
            f.write('%10d %12.5f \n' % (key, value))


def load_dict_orcaS(fileName: Path) -> dict:
    IsExist(fileName)
    lines: list = open(fileName, "r").readlines()
    Data: dict = {}
    for x in lines:
        Data[int(x.split()[0])] = float(x.split()[1])
    return Data


def IsExist(fileName: Path) -> None:
    IsExists: bool = Path(fileName).exists()
    if not IsExists:
        print(f"{fileName} the file is not exist ...")
        print("    Exit and close the program !!! ")
        ic()
        raise FileNotFoundError(f"{fileName} The file is not Exist ...")


def IsExist_return_bool(fileName: Path) -> bool:
    IsExists: bool = Path(fileName).exists()
    if IsExists:
        return True
    else:
        print("="*80)
        print(f"{fileName} the file is not exist ...")
        print("="*80)
        return False


def program_IsExist(ProgramName: str) -> bool:
    from shutil import which
    if which(ProgramName):
        return True
    else:
        print(ProgramName, " the program is not exist ...")
        print(" Exit and close the program !!! ")
        ic()
        raise ValueError(" the program is not Exist ...")


def save_figure(fileName="nmrplot") -> None:
    import matplotlib.pyplot as plt
    plt.savefig(fileName + ".pdf", dpi=300)
    plt.savefig(fileName + ".svg")
