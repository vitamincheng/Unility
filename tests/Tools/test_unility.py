from pathlib import Path
import pytest


def test_unility_IsExists_DirFileName_miss_args():
    from censo_ext.Tools.utility import IsExists_DirFileName

    FileName = Path("tests/data/34.Ergocalciferol/04.Hydrogen/test.dat")
    with pytest.raises(FileNotFoundError) as e:
        a, b = IsExists_DirFileName(FileName)
    assert str(e.value) == f"  The file {FileName} is not Exist ..."


def test_unility_IsExists_DirFileName():
    from censo_ext.Tools.utility import IsExists_DirFileName
    assert IsExists_DirFileName(
        Path("tests/data/34.Ergocalciferol/04.Hydrogen/.anmrrc")) == (Path("tests/data/34.Ergocalciferol/04.Hydrogen"), ".anmrrc")
    a, b = IsExists_DirFileName(
        Path("tests/data/34.Ergocalciferol/04.Hydrogen/anmr.dat"))
    assert a == Path("tests/data/34.Ergocalciferol/04.Hydrogen")
    assert b == "anmr.dat"


@pytest.mark.parametrize(argnames="input_bool,input_str",
                         argvalues=[(True, "123.5678"), (True, " 123.5678"),
                                    (True, "123.5678 "), (True, " 123.5678 "),
                                    (False, " test "), (False, " 123.56EA ")])
def test_unility_float(input_bool: bool, input_str: str):
    from censo_ext.Tools.utility import function_is_float
    assert function_is_float(input_str) == input_bool


@pytest.mark.parametrize(argnames="input_bool,input_str",
                         argvalues=[(True, "456"), (True, " 123 "),
                                    (False, "456.1234"), (False, " test ")])
def test_unility_int(input_bool: bool, input_str: str):
    from censo_ext.Tools.utility import function_is_int
    assert function_is_int(input_str) == input_bool


@pytest.mark.parametrize(argnames="input_bool,input_Path",
                         argvalues=[(True, Path("tests/Tools/test_unility.py")),
                                    (False, Path("tests/Tools/test_unility1.py"))])
def test_unility_IsExistReturnBool(input_bool: bool, input_Path: Path):
    from censo_ext.Tools.utility import IsExist_return_bool
    assert IsExist_return_bool(input_Path) == input_bool


def test_unility_ProgramIsExist():
    from censo_ext.Tools.utility import program_IsExist
    assert program_IsExist("xtb")
    assert program_IsExist("orca")
    with pytest.raises(ValueError) as e:
        program_IsExist("kkk")
    assert str(e.value) == " the program is not Exist ..."


def test_unility_unilityIsExist():
    from censo_ext.Tools.utility import IsExist
    assert IsExist(Path("tests/data/crest_conformers.xyz")) is None
    fileName: str = "kkk.xyz"
    with pytest.raises(FileNotFoundError) as e:
        IsExist(Path(fileName))
    assert str(e.value) == f"  The file {fileName} is not Exist ..."


def test_unility_move_file():
    from censo_ext.Tools.utility import move_file
    import os
    source = Path("/tmp/test_source.txt")
    destination = Path("/tmp/test_destination.txt")

    with open(source, 'w') as f:
        f.write("Hello World!")

    move_file(source, destination)
    assert os.path.exists(destination)


def test_unility_copy_file():
    from censo_ext.Tools.utility import copy_file
    import os
    source = Path("/tmp/test_source.txt")
    destination = Path("/tmp/test_destination.txt")

    with open(source, 'w') as f:
        f.write("Hello World!")

    copy_file(source, destination)
    assert os.path.exists(destination)


def test_unility_delete_all_files():
    from censo_ext.Tools.utility import delete_all_files
    import os
    files_to_delete = [Path("/tmp/file1.txt"), Path("/tmp/file2.txt")]

    for file in files_to_delete:
        with open(file, 'w') as f:
            f.write("Hello World!")

    delete_all_files(*files_to_delete)
    for file in files_to_delete:
        assert not os.path.exists(file)


def test_unility_delete_file_bool():
    from censo_ext.Tools.utility import delete_file_bool
    import os
    file = Path("/tmp/test_file.txt")

    with open(file, 'w') as f:
        f.write("Hello World!")

    assert delete_file_bool(file) is True

    assert not os.path.exists(file)


def test_unlity_save_dict():
    from censo_ext.Tools.utility import save_dict
    data = {1: 2.3456789, 2: 3.4567890}

    save_dict(Path("/tmp/test_dict.txt"), data)

    with open(Path("/tmp/test_dict.txt")) as f:
        lines = f.readlines()

        assert len(lines) == 2
        assert lines[0].strip() == "1.00000  2.34568e+00"
        assert lines[1].strip() == "2.00000  3.45679e+00"


def test_unility_save_dict_orcaS():
    from censo_ext.Tools.utility import save_dict_orcaS
    data = {1: 2.3456789}

    save_dict_orcaS(Path("/tmp/test_dict.txt"), data)

    with open(Path("/tmp/test_dict.txt")) as f:
        lines = f.readlines()

        assert len(lines) == 1
        assert lines[0].strip() == "1      2.34568"


def test_unility_load_dict_orcaS():
    from censo_ext.Tools.utility import load_dict_orcaS
    data = {1: 2.3456789}

    with open(Path("/tmp/test_dict.txt"), 'w') as f:
        f.write('1      2.3456789\n')

    loaded_data = load_dict_orcaS(Path("/tmp/test_dict.txt"))

    assert loaded_data == data
