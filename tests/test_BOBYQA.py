
#!/usr/bin/env python3
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
RefDat: Path = Path("1r_h.dat")


def test_cregen_miss_args():
    args_x: dict = {}
    args = argparse.Namespace(**args_x)
    with pytest.raises(SystemExit) as e:
        BOBYQA.main(args)
    assert e.type == SystemExit
    assert e.value.code == 2  # for argparse error


def test_BOBYQA_init():
    args_x: dict = {"auto": True, "average": False, "dir": DirName,
                    "bobyqa": False, "mf": 500, "lw": None, "thr": None, "json": None, "thrab": 0.025,
                    "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
    args = argparse.Namespace(**args_x)
    anmr.main(args)

    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    # Create orcaS_BOBYQA file
    BOBYQA.main(args2)


def test_BOBYQA_block_file():
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    # run block orcaS_BOBYQA file
    BOBYQA.main(args2)


def test_BOBYQA_single():
    import shutil
    shutil.copyfile(DirName / Path("orcaS-BOBYQA.out"),
                    DirName / Path("Average/NMR/orcaS-BOBYQA.out"))

    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": None}
    args2 = argparse.Namespace(**args_y)
    BOBYQA.main(args2)


@pytest.mark.skipif(_system == "Darwin", reason="crest only work under linux")
def test_BOBYQA_single_external_prog():
    args_y: dict = {"dir": DirName, "ref": RefDat, "limit": 0.20, "prog": True}
    args2 = argparse.Namespace(**args_y)
    BOBYQA.main(args2)

    # compare = ""
    # if _system == "Linux":
    #    compare = f"tests/compare/cregen_cluster.xyz"
    # elif _system == "Darwin":  # No crest under Darwin system
    #    compare = f"tests/compare/cregen_cluster_Darwin.xyz"
    # else:
    #    pytest.raises(
    #        ValueError, match="OS system only can run under Darwin or Linux")
    # assert filecmp.cmp(args.out, compare) == True
    # os.remove(args.out)
    # os.remove("isomers.out")


if __name__ == "__main__":
    import cProfile
    cProfile.run("test_BOBYQA_init()", sort="cumtime")
    cProfile.run("test_BOBYQA_single()", sort="cumtime")
