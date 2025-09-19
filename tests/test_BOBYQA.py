#!/usr/bin/env python
import pytest
from pathlib import Path
import argparse
import censo_ext.BOBYQA as BOBYQA
import censo_ext.anmr as anmr
import filecmp
import platform

_system = platform.system()
DirName: Path = Path(
    "tests/data/31.Cyclohexanone/03.Censo_for_Hydrogen_(revTPSS)")
DirCompare: Path = Path("tests/compare/BOBYQA")
RefDat: Path = Path("1r_h.dat")
args_Normal: dict = {"dir": DirName, "ref": RefDat, "mf": 500, "lw": 1,
                     "limit": 0.20, "prog": None, "verbose": False}


def anmr_init():
    # Create Average/NMR/orcaS.out and Average/NMR/orcaS_BOBYQA.out file
    args_x: dict = {"auto": True, "average": False, "dir": DirName, "bobyqa": False, "mf": 500,
                    "lw": None, "thr": None, "json": None, "thrab": 0.025, "tb": 4, "mss": 9, "verbose": False,
                    "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    args = argparse.Namespace(**args_x)
    anmr.main(args)


def BOBYQA_init():
    args2 = argparse.Namespace(**args_Normal)

    with pytest.raises(SystemExit) as e:
        BOBYQA.main(args2)
    assert e.type is SystemExit
    assert e.value.code == 0  # for normal exit(0)


def BOBYQA_final_remove_files():
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(DirName / "output.dat",
                     DirName / "peaks.json",
                     DirName / "Average/NMR/orcaA.out",
                     DirName / "Average/NMR/orcaJ.out",
                     DirName / "Average/NMR/orcaS.out",
                     DirName / "Average/NMR/orcaS-BOBYQA.out",)


def test_BOBYQA_miss_all_args():
    args_x: dict = {}
    with pytest.raises(SystemExit) as e:
        anmr.main(argparse.Namespace(**args_x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error

    args_x: dict = {}
    with pytest.raises(SystemExit) as e:
        BOBYQA.main(argparse.Namespace(**args_x))
    assert e.type is SystemExit
    assert e.value.code == 2  # for argparse error


def test_BOBYQA_miss_orcaS_file():
    #
    # not run anmr_BOBYQA_init() and directly run BOBYQA.py
    BOBYQA_final_remove_files()
    with pytest.raises(FileNotFoundError) as e:
        BOBYQA.main(argparse.Namespace(**args_Normal))
    assert str(e.value) == str(
        DirName/Path("Average/NMR/orcaS.out")) + " is not exist !!!"


def test_BOBYQA_blank_file():
    # run block orcaS_BOBYQA file
    anmr_init()
    BOBYQA_init()


def test_BOBYQA_blank_file_final_remove_files():
    BOBYQA_final_remove_files()


@pytest.mark.slow
def test_BOBYQA_single():
    anmr_init()
    BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))
    BOBYQA.main(argparse.Namespace(**args_Normal))
    assert filecmp.cmp(DirName / Path("output.dat"),
                       DirCompare / Path("orcaS-BOBYQA-anmrpy.dat"))
    BOBYQA_final_remove_files()


@pytest.mark.slow
@pytest.mark.skipif(_system == "Darwin", reason="anmr only work under linux CLI")
def test_BOBYQA_single_external_prog(monkeypatch):
    anmr_init()
    BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))

    from io import StringIO
    monkeypatch.setattr('sys.stdin', StringIO('Y\n'))
    BOBYQA.main(argparse.Namespace(**args_Normal))
    assert filecmp.cmp(DirName / Path("anmr.dat"),
                       DirCompare / Path("orcaS-BOBYQA-anmr.dat"))
    BOBYQA_final_remove_files()
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(DirName / "anmr.dat",
                     DirName / "anmr.out",
                     DirName / "tmpanmr.1",
                     DirName / "tmpanmr_frag.av",
                     DirName / "tmpanmr_full.av")


@pytest.mark.slow
def test_BOBYQA_group():
    anmr_init()
    BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA-group.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))
    BOBYQA.main(argparse.Namespace(**args_Normal))
    assert filecmp.cmp(DirName / Path("output.dat"),
                       DirCompare / Path("orcaS-BOBYQA-group-anmrpy.dat"))
    BOBYQA_final_remove_files()


@pytest.mark.slow
@pytest.mark.skipif(_system == "Darwin", reason="anmr only work under linux CLI")
def test_BOBYQA_group_external_prog(monkeypatch):
    anmr_init()
    BOBYQA_init()
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA-group.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))
    from io import StringIO
    monkeypatch.setattr('sys.stdin', StringIO('Y\n'))
    BOBYQA.main(argparse.Namespace(**args_Normal))
    assert filecmp.cmp(DirName / Path("anmr.dat"),
                       DirCompare / Path("orcaS-BOBYQA-group-anmr.dat"))
    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(DirName / "anmr.dat", DirName / "anmr.out")


if __name__ == "__main__":
    import cProfile
    from pstats import Stats

    pr = cProfile.Profile()
    pr.enable()

    test_BOBYQA_single()

    pr.disable()
    stats = Stats(pr)
    stats.sort_stats('cumtime').print_stats(50)

# if __name__ == "__main__":
# import cProfile
# cProfile.run("test_BOBYQA_init()", sort="cumtime")
# cProfile.run("test_BOBYQA_single()", sort="cumtime")


# For all test_fuction
#   pytest -v -s tests -m "not slow"
#   pytest -v -s tests -m slow

# For slow test_function
#  (--pdb is debug,-v verbose, -s display print)
#   pytest -v -s --pdb tests/test_BOBYQA.py::test_BOBYQA_single_prog
#   pytest -v -s --pdb tests/test_BOBYQA.py::test_BOBYQA_single_external_prog
#   pytest -v -s --pdb tests/test_BOBYQA.py::test_BOBYQA_group_prog
#   pytest -v -s --pdb tests/test_BOBYQA.py::test_BOBYQA_group_external_prog
