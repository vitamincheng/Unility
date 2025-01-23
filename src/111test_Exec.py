#! /usr/bin/env python3
from numpy import ndarray

if __name__ == "__main__":

    import argparse
    from icecream import ic
    content: dict[int, str] = {
        1:   "anmr.py",
        2:   "BOBYQA.py",
        3:   "cregen.py",
        4:   "datNormailzed.py",
        101: "ensoAnalysis.py",
        102: "ensoGenFlexible.py",
        5:   "FactorAnalysis.py",
        6:   "FactorCompare.py",
        7:   "molclus_orca.py",
        8:   "molclus_xtb.py",
        9:   "molManipulate.py",
        10:  "nmrplotSJ.py",
        103: "plot_1D_DEPT.py",
        11:  "plot_1D.py",
        12:  "plot_2D_exp.py",
        13:  "plot_2D_sim.py",
        104: "xyzGenFlexible.py",
        14:  "xyzReturnOandZ.py",
        15:  "xyzSerial.py",
        16:  "xyzSplit.py",
        17:  "xyzTranslate.py",
        99:  "quit", }

    ic(content)
    nCH: int = int(input())
    while nCH != 99:

        print("="*80)
        print(content[nCH])
        print("="*80)

        ################################################
        # anmr.py
        ################################################
        if nCH == 1:
            import anmr
            import argparse
            x: dict = {"auto": True, "average": True, "dir": "../tests/04.Hydrogen",
                       "bobyqa": True, "mf": 500, "thr": None, "json": [-1], "thrab": 0.025,
                       "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
            peaks: ndarray = anmr.main(argparse.Namespace(**x))
            ic(peaks.shape)

            x: dict = {"auto": True, "average": True, "dir": "../tests/07.Carbon",
                       "bobyqa": True, "mf": 500, "thr": None, "json": [-1], "thrab": 0.025,
                       "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
            peaks: ndarray = anmr.main(argparse.Namespace(**x))
            ic(peaks.shape)

            x: dict = {"auto": True, "average": True, "dir": "../tests/04.Hydrogen",
                       "bobyqa": False, "mf": 500, "thr": None, "json": None, "thrab": 0.025,
                       "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
            anmr.main(argparse.Namespace(**x))

            x: dict = {"auto": True, "average": True, "dir": "../tests/07.Carbon",
                       "bobyqa": False, "mf": 500, "thr": None, "json": None, "thrab": 0.025,
                       "tb": 4, "mss": 9, "cutoff": 0.001, "show": False, "start": None, "end": None, "out": "output.dat"}
            anmr.main(argparse.Namespace(**x))
            # python3 anmr.py -d ../tests/04.Hydrogen --auto -av -j -1
            # python3 anmr.py -d ../tests/07.Carbon --auto -av -j -1
            # python3 anmr.py -d ../tests/04.Hydrogen --auto -av
            # python3 anmr.py -d ../tests/07.Carbon --auto -av

        ################################################
        # cregen.py
        ################################################
        if nCH == 3:
            import cregen
            import argparse
            x: dict = {"file": "../tests/crest_conformers.xyz", "rthr": 0.175,
                       "bthr": 0.03, "ethr": 0.15, "ewin": 4, "out": "cluster.xyz"}
            cregen.main(argparse.Namespace(**x))

        ################################################
        # datNormalized.py
        ################################################
        if nCH == 4:
            import datNormalized
            import argparse

            x = {"file": "../tests/04.Hydrogen/output.dat",
                 "start": -5, "end": 15, "dpi": 10000, "out": "output.dat"}
            datNormalized.main(argparse.Namespace(**x))

            x = {"file": "../tests/07.Carbon/output.dat",
                 "start": -5, "end": 15, "dpi": 10000, "out": "output.dat"}
            datNormalized.main(argparse.Namespace(**x))
            # python3 datNormalized.py -i ../tests/04.Hydrogen/output.dat
            # python3 datNormalized.py -i ../tests/07.Carbon/output.dat

        ################################################
        # FactorAnalysis.py
        ################################################
        if nCH == 5:
            import FactorAnalysis
            import argparse
            x = {"auto": True, "file": "../tests/crest_conformers.xyz",
                 "Analysis": True, "factor": None, "opt": None, "Filter": False}
            FactorAnalysis.main(args=argparse.Namespace(**x))
            x = {"auto": True, "file": "../tests/crest_conformers.xyz",
                 "Analysis": True, "factor": None, "opt": True, "Filter": False}
            FactorAnalysis.main(args=argparse.Namespace(**x))
            x = {"auto": True, "file": "../tests/crest_conformers.xyz", "thr": 2, "remove_idx": None, "add_idx": None,
                 "bond_broken": [52, 55], "ignore_Hydrogen": True, "Analysis": False, "factor": None, "opt": None, "Filter": True}
            FactorAnalysis.main(argparse.Namespace(**x))
            # python3 FactorAnalysis.py -A -i ../tests/crest_conformers.xyz --opt
            # python3 FactorAnalysis.py -F -i ../tests/crest_conformers.xyz -bb 52 55 -nh

        ################################################
        # molclus_orca.py
        ################################################
        if nCH == 7:
            import molclus_orca
            import argparse
            x = {"file": "../tests/crest_conformers.xyz", "template": "template.inp", "remove": True,
                 "chrg": 0, "uhf": 1, "out": "isomers.xyz"}
            molclus_orca.main(argparse.Namespace(**x))
            x = {"file": "../tests/crest_conformers.xyz", "template": "template.inp", "remove": True,
                 "chrg": 0, "uhf": 1, "out": "isomers.xyz"}
            molclus_orca.main(argparse.Namespace(**x))

            # python3 molclus_xtb.py -i ../tests/crest_conformers.xyz
            # python3 molclus_xtb.py -i ../tests/crest_conformers.xyz --opt

        ################################################
        # molclus_xtb.py
        ################################################
        if nCH == 8:
            import molclus_xtb
            import argparse
            x = {"file": "../tests/crest_conformers.xyz", "method": "gfn2",
                 "chrg": 0, "uhf": 1, "out": "isomers.xyz", "alpb": None, "gbsa": None, "opt": False}
            molclus_xtb.main(argparse.Namespace(**x))
            x = {"file": "../tests/crest_conformers.xyz", "method": "gfn2",
                 "chrg": 0, "uhf": 1, "out": "isomers.xyz", "alpb": None, "gbsa": None, "opt": True}
            molclus_xtb.main(argparse.Namespace(**x))

            # python3 molclus_xtb.py -i ../tests/crest_conformers.xyz
            # python3 molclus_xtb.py -i ../tests/crest_conformers.xyz --opt

        ################################################
        # molManipulate.py
        ################################################
        if nCH == 9:
            import molManipulate
            import argparse
            x = {"separate": "../tests/crest_conformers3.xyz"}
            molManipulate.main(argparse.Namespace(**x))

            x = {"separate": None, "out": "output.xyz", "merge": [
                "Separation/1.xyz", "Separation/2.xyz"]}
            molManipulate.main(argparse.Namespace(**x))
            # python3 molManipulate.py -s ../tests/crest_conformers3.xyz
            # python3 molManipulate.py -m Separation/1.xyz Separation/2.xyz

        ################################################
        # xyzReturnOandZ.py
        ################################################
        if nCH == 14:
            import xyzReturnOandZ
            import argparse
            x: dict = {"file": "../tests/crest_conformers.xyz",
                       "atom": [30, 45, 47], "print": False, "replace": False, "out": "output.xyz"}
            xyzReturnOandZ.main(argparse.Namespace(**x))
            x: dict = {"file": "../tests/crest_conformers.xyz", "auto": True,
                       "atom": None, "print": False, "replace": False, "out": "output.xyz"}
            xyzReturnOandZ.main(argparse.Namespace(**x))
            # python3 xyzReturnOandZ.py -i Tools/crest_conformers.xyz -a 30 45 47
            # python3 xyzReturnOandZ.py -i Tools/crest_conformers.xyz --auto

        ################################################
        # xyzSerial.py
        ################################################
        if nCH == 15:
            import xyzSerial
            import argparse
            x = {"file": "../tests/crest_conformers.xyz", "new": True,
                 "keep": False, "print": False, "out": "output.xyz"}
            xyzSerial.main(argparse.Namespace(**x))
            x = {"file": "../tests/crest_conformers.xyz", "new": False,
                 "keep": True, "print": False, "out": "output.xyz"}
            xyzSerial.main(argparse.Namespace(**x))
            x = {"file": "../tests/crest_conformers.xyz", "new": True,
                 "keep": False, "print": True, "out": "output.xyz"}
            xyzSerial.main(argparse.Namespace(**x))
            # python3 xyzSerial.py -i Tools/crest_conformers.xyz -n
            # python3 xyzSerial.py -i Tools/crest_conformers.xyz -k
            # python3 xyzSerial.py -i Tools/crest_conformers.xyz -n -p

        ################################################
        # xyzSplit.py
        ################################################
        if nCH == 16:
            import xyzSplit
            import argparse
            x = {"file": "../tests/crest_conformers1.xyz",
                 "atom": [52, 55], "cut": 12, "print": False, "out": "output.xyz"}
            xyzSplit.main(argparse.Namespace(**x))
            x = {"file": "../tests/crest_conformers1.xyz",
                 "atom": [52, 55], "cut": 12, "print": True, "out": "output.xyz"}
            xyzSplit.main(argparse.Namespace(**x))
            #   python3 xyzSplit.py -i Tools/crest_conformers.xyz -a 52 55 -c 3
            #   python3 xyzSplit.py -i Tools/crest_conformers.xyz -a 52 55 -c 3 -p

        ################################################
        # xyzTranslate.py
        ################################################
        if nCH == 17:
            import xyzTranslate
            import argparse
            x = {"file": "../tests/crest_conformers.xyz",
                 "move": [5, 0, 0], "cut": None, "out": "output.xyz"}
            xyzTranslate.main(argparse.Namespace(**x))
            x = {"file": "../tests/crest_conformers1.xyz",
                 "move": [5, 0, 0], "cut": 10, "out": "output.xyz"}
            xyzTranslate.main(argparse.Namespace(**x))
            x = {"file": "output.xyz", "move": [
                0, 0, 5], "cut": 10, "out": "output2.xyz"}
            xyzTranslate.main(argparse.Namespace(**x))
            # python3 xyzTranslate.py -i Tools/crest_conformers.xyz -m 5 0 0
            # python3 xyzTranslate.py -i Tools/crest_conformers.xyz -m 5 0 0 -c 10
            # python3 xyzTranslate.py -i output.xyz -m 0 0 5 -c 10 -o output2.xyz

        ic(content)
        nCH = int(input())
