#! /usr/bin/env python3
import numpy as np
from icecream import ic


def main():
    print("test")
    cross = np.genfromtxt("plot_2D_exp.peaks", names=True)
    ic(cross)
    Carbon = np.genfromtxt("plot_1D_DEPT.peaks", names=True)
    ic(Carbon)
    Hydrogen = np.genfromtxt("plot_1D.peaks", names=True)
    ic(Hydrogen)
    ic(Hydrogen['Area'])
    one_Area = Hydrogen['Area'].sum()/43

    ic(Hydrogen['Area']/one_Area)


if __name__ == "__main__":
    main()
