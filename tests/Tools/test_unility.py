from pathlib import Path
import pytest


def test_Unility_IsExistsDirFileName():
    from censo_ext.Tools.utility import IsExistsDirFileName
    a, b = IsExistsDirFileName(Path("tests/data/04.Hydrogen/anmr.dat"))
    assert a == Path("tests/data/04.Hydrogen")
    assert b == "anmr.dat"


def test_unility_file():
    from censo_ext.Tools.utility import move_file, copy_file, delete_all_files, delete_file_bool
    copy_file(Path("tests/data/crest_conformers.xyz"),
              Path("tests/data/temp.xyz"))
    move_file(Path("tests/data/temp.xyz"), Path("tests/data/temp1.xyz"))
    copy_file(Path("tests/data/temp1.xyz"), Path("tests/data/temp2.xyz"))
    delete_all_files(Path("tests/data/temp1.xyz"),
                     Path("tests/data/temp2.xyz"))
    copy_file(Path("tests/data/crest_conformers.xyz"),
              Path("tests/data/temp.xyz"))
    assert False == delete_file_bool(Path("tests/data/temp1.xyz"))
    assert True == delete_file_bool(Path("tests/data/temp.xyz"))


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
    from censo_ext.Tools.utility import is_exist_return_bool
    from pathlib import PosixPath
    assert is_exist_return_bool(Path("tests/Tools/test_unility.py")) == True
    assert is_exist_return_bool(
        Path("tests/Tools/test_unility1.py")) == False


def test_unility_IsExistsDirFileName():
    from censo_ext.Tools.utility import IsExistsDirFileName
    assert IsExistsDirFileName(
        Path("tests/data/04.Hydrogen/.anmrrc")) == (Path("tests/data/04.Hydrogen"), ".anmrrc")


def test_unility_ProgramIsExist():
    from censo_ext.Tools.utility import program_is_exist
    assert program_is_exist("xtb") == True
    assert program_is_exist("orca") == True
    with pytest.raises(SystemExit) as e:
        program_is_exist("kkk")
    assert e.type == SystemExit
    assert e.value.code == 1


def test_unility_is_exist():
    from censo_ext.Tools.utility import is_exist
    assert is_exist(Path("tests/data/crest_conformers.xyz")) == None
    with pytest.raises(SystemExit) as e:
        is_exist(Path("kkk.xyz"))
    assert e.type == SystemExit
    assert e.value.code == 1

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
