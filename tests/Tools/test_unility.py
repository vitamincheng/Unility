import pathlib
import argparse
from pathlib import Path
import pytest


def test_unility_float():
    from censo_ext.Tools.utility import function_is_float
    assert function_is_float("123.5678") == True
    assert function_is_float(" 123.5678") == True
    assert function_is_float("123.5678 ") == True
    assert function_is_float(" test ") == False


def test_unility_int():
    from censo_ext.Tools.utility import function_is_int
    assert function_is_int("123.5678") == False
    assert function_is_int(" 123 ") == True
    assert function_is_int(" test ") == False


def test_unility_IsExistReturnBool():
    from censo_ext.Tools.utility import IsExistReturnBool
    from pathlib import PosixPath
    assert IsExistReturnBool(Path("tests/Tools/test_unility.py")) == True
    assert IsExistReturnBool(
        Path("tests/Tools/test_unility1.py")) == False


def test_unility_IsExistsDirFileName():
    from censo_ext.Tools.utility import IsExistsDirFileName
    assert IsExistsDirFileName(
        Path("tests/data/04.Hydrogen/.anmrrc")) == (Path("tests/data/04.Hydrogen"), ".anmrrc")


def test_unility_ProgramIsExist():
    from censo_ext.Tools.utility import ProgramIsExist
    assert ProgramIsExist("xtb") == True
    assert ProgramIsExist("orca") == True


#
#
#
#
#
# def example():
#    print("foo")
#
#
# def test_spam(capsys):
#    example()
#    captured = capsys.readouterr()
#    assert captured.out == "foo\n"
