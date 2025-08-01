#!/usr/bin/env python3
import argparse
__doc__ = """
             [07.14.2023]   vitamin.cheng@gmail.com
                                            
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation.

For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

from icecream import ic
import copy
import numpy as np
import numpy.typing as npt
from pathlib import Path
import censo_ext.Tools.Parameter as Parameter
from censo_ext.Tools.xyzfile import GeometryXYZs

# METHOD_KABSCH = "kabsch"
# ROTATION_METHODS = [METHOD_KABSCH]

NAMES_ELEMENT = {value: key for key, value in Parameter.ELEMENT_NAMES.items()}


# class RmsdCallable(Protocol):
#    def __call__(self, P: ndarray, Q: ndarray, **kwargs: Any,) -> float | None:
#        """
#        Protocol for a rotation callable function
#
#        return:
#            RMSD after rotation
#        """
#

def str_atom(atom: int) -> str:
    """
    Convert atom type from integer to string
    """
    return Parameter.ELEMENT_NAMES[atom]


def int_atom(atom: str) -> int:
    """
    Convert atom type from string to integer
    """
    atom = atom.capitalize().strip()
    return NAMES_ELEMENT[atom]


def rmsd(P: npt.NDArray[np.float64], Q: npt.NDArray[np.float64], idx_atom: list[int], **kwargs) -> tuple[dict[int, float], float]:
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    diff: npt.NDArray[np.float64] = P - Q
    atom_square: dict[int, float] = {}
    coord_square_total: float = 0
    for idx0, x in enumerate(idx_atom):
        coord_square: float = float((diff[idx0]*diff[idx0]).sum())
        # ic(coord_square)
        if __name__ == "__main__":
            print(f"{x:>5}", end=" ")
            print(f"{coord_square:>10.5f}")
        atom_square[x] = coord_square
        coord_square_total += coord_square
    return atom_square, float(np.sqrt(coord_square_total / P.shape[0]))


def kabsch_rmsd(P: npt.NDArray, Q: npt.NDArray, idx_atom1: list,
                translate: bool = False) -> tuple[dict[int, float], float]:
    """
    Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.
    An optional vector of weights W may be provided.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    translate : bool
        Use centroids to translate vector P and Q unto each other.

    Returns
    -------
    rmsd : float
        root-mean squared deviation
    """

    if translate:
        Q = Q - centroid(Q)
        P = P - centroid(P)

    P = kabsch_rotate(P, Q)
    A, B = rmsd(P, Q, idx_atom1)
    return A, B


def kabsch_rotate(P: npt.NDArray[np.float64], Q: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    P : array
        (N,D) matrix, where N is points and D is dimension,
        rotated

    """
    U: npt.NDArray[np.float64] = kabsch(P, Q)

    # Rotate P
    Res: npt.NDArray[np.float64] = np.dot(P, U)
    return Res


def kabsch(P: npt.NDArray[np.float64], Q: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    The algorithm works in three steps:
    - a centroid translation of P and Q (assumed done before this function
      call)
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    U : matrix
        Rotation matrix (D,D)
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U: npt.NDArray[np.float64] = np.dot(V, W)

    return U


# def kabsch_weighted(P: ndarray, Q: ndarray, W: Optional[ndarray] = None) -> Tuple[ndarray, ndarray, float]:
#    """
#    Using the Kabsch algorithm with two sets of paired point P and Q.
#    Each vector set is represented as an NxD matrix, where D is the
#    dimension of the space.
#    An optional vector of weights W may be provided.
#
#    Note that this algorithm does not require that P and Q have already
#    been overlayed by a centroid translation.
#
#    The function returns the rotation matrix U, translation vector V,
#    and RMS deviation between Q and P', where P' is:
#
#        P' = P * U + V
#
#    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm
#
#    Parameters
#    ----------
#    P : array
#        (N,D) matrix, where N is points and D is dimension.
#    Q : array
#        (N,D) matrix, where N is points and D is dimension.
#    W : array or None
#        (N) vector, where N is points.
#
#    Returns
#    -------
#    U    : matrix
#           Rotation matrix (D,D)
#    V    : vector
#           Translation vector (D)
#    RMSD : float
#           Root mean squared deviation between P and Q
#    """
#    # Computation of the weighted covariance matrix
#    CMP = np.zeros(3)
#    CMQ = np.zeros(3)
#    C = np.zeros((3, 3))
#    if W is None:
#        W = np.ones(len(P)) / len(P)
#    W = np.array([W, W, W]).T
#    # NOTE UNUSED psq = 0.0
#    # NOTE UNUSED qsq = 0.0
#    iw = 3.0 / W.sum()
#    n = len(P)
#    for i in range(3):
#        for j in range(n):
#            for k in range(3):
#                C[i, k] += P[j, i] * Q[j, k] * W[j, i]
#    CMP = (P * W).sum(axis=0)
#    CMQ = (Q * W).sum(axis=0)
#    PSQ = (P * P * W).sum() - (CMP * CMP).sum() * iw
#    QSQ = (Q * Q * W).sum() - (CMQ * CMQ).sum() * iw
#    C = (C - np.outer(CMP, CMQ) * iw) * iw
#
#    # Computation of the optimal rotation matrix
#    # This can be done using singular value decomposition (SVD)
#    # Getting the sign of the det(V)*(W) to decide
#    # whether we need to correct our rotation matrix to ensure a
#    # right-handed coordinate system.
#    # And finally calculating the optimal rotation matrix U
#    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
#    V, S, W = np.linalg.svd(C)
#    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0  # type: ignore
#
#    if d:
#        S[-1] = -S[-1]
#        V[:, -1] = -V[:, -1]
#
#    # Create Rotation matrix U, translation vector V, and calculate RMSD:
#    U = np.dot(V, W)  # type: ignore
#    msd = (PSQ + QSQ) * iw - 2.0 * S.sum()
#    if msd < 0.0:
#        msd = 0.0
#    rmsd_ = np.sqrt(msd)
#    V = np.zeros(3)
#    for i in range(3):
#        t = (U[i, :] * CMQ).sum()
#        V[i] = CMP[i] - t
#    V = V * iw
#    return U, V, rmsd_
#
#
# def kabsch_weighted_rmsd(P: ndarray, Q: ndarray, W: Optional[ndarray] = None) -> float:
#    """
#    Calculate the RMSD between P and Q with optional weighhts W
#
#    Parameters
#    ----------
#    P : array
#        (N,D) matrix, where N is points and D is dimension.
#    Q : array
#        (N,D) matrix, where N is points and D is dimension.
#    W : vector
#        (N) vector, where N is points
#
#    Returns
#    -------
#    RMSD : float
#    """
#    _, _, w_rmsd = kabsch_weighted(P, Q, W)
#    return w_rmsd


def centroid(X: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Centroid is the mean position of all the points in all of the coordinate
    directions, from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid

    C = sum(X)/len(X)

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : ndarray
        centroid
    """
    return X.mean(axis=0)


def get_Coordinates(xyzfile, idx0) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    '''
    Read xyz file to data 
    idx is Serial No. (from 0) in xyz file
    '''
    Names: dict[int, str] = xyzfile.Sts[idx0].names
    element: npt.NDArray[np.int64] = np.array(
        [int_atom(atom) for atom in Names.values()])
    V: npt.NDArray[np.float64] = np.array(
        xyzfile.Sts[idx0].coord, dtype=np.float64)
    return element, V


def cal_RMSD_xyz(xyzfile: GeometryXYZs, idx_p: int, idx_q: int, args: argparse.Namespace) -> tuple[dict[int, float], float]:
    '''
    Read xyz file and calculate rmsd
    xyzfile is class ClassGeometryXYZs
    idx_p is from 1 to Serial No.
    idx_q is from 1 to Serial No.
    '''
    xyz_tmp: Path = Path(".tmp.xyz")
    idx_p -= 1
    idx_q -= 1
    p_all_atoms: npt.NDArray[np.int64]
    q_all_atoms: npt.NDArray[np.int64]
    p_all: npt.NDArray[np.float64]
    q_all: npt.NDArray[np.float64]
    p_all_atoms, p_all = get_Coordinates(xyzfile, idx_p)
    q_all_atoms, q_all = get_Coordinates(xyzfile, idx_q)

    idx_atom1: npt.NDArray[np.int64] = np.array([], dtype=np.int64)

    if p_all.shape[0] != q_all.shape[0]:
        print("error: Structures not same size")
        ic()
        raise ValueError("error: Structures not same size")

    # Typing
    from typing import Union
    index: Union[set[int], list[int], npt.NDArray[np.int64]]

    # Set local view
    p_view: None | npt.NDArray[np.int64] = None
    q_view: None | npt.NDArray[np.int64] = None

    if args.ignore_Hydrogen:
        assert type(p_all_atoms[0]) is not str
        assert type(q_all_atoms[0]) is not str
        for idx in range(len(p_all_atoms)):
            if p_all_atoms[idx] != 1:
                idx_atom1 = np.append(idx_atom1, [idx+1])

        p_view = np.where(p_all_atoms != 1)  # type: ignore
        q_view = np.where(q_all_atoms != 1)  # type: ignore

    if args.bond_broken:
        if args.ignore_Hydrogen:
            xyzfile.set_filename(xyz_tmp)
            xyzfile.method_save_xyz([idx_p])
            args_x: dict = {"file": xyz_tmp, "bond_broken": [
                args.bond_broken[0], args.bond_broken[1]], "print": False, "debug": False}
            from censo_ext.Tools.topo import Topo
            Sts_topo: Topo = Topo(args_x["file"])
            idx_atom1 = np.array(Sts_topo.method_broken_bond(
                argparse.Namespace(**args_x)))
            index = idx_atom1-1
            p_view, q_view = index, index

        else:
            print("Only support under ignore Hydrogen condition ")
            ic()
            raise ValueError("Only support under ignore Hydrogen condition ")
    else:
        pass

    if args.remove_idx:
        if args.ignore_Hydrogen:
            pass
        else:
            idx_atom1 = np.arange(len(p_all_atoms), dtype=np.int64)+1

        idx_atom1 = np.setdiff1d(idx_atom1, args.remove_idx)

        args.remove_idx = np.array(args.remove_idx)-1
        index = idx_atom1-1

        p_view, q_view = index, index

    elif args.add_idx:
        if args.ignore_Hydrogen:
            idx_atom1 = np.union1d(idx_atom1, args.add_idx)
        else:
            idx_atom1 = args.add_idx
        args.add_idx = idx_atom1-1

        p_view, q_view = args.add_idx, args.add_idx

    # Set local view
    if p_view is None:
        p_coord: npt.NDArray[np.float64] = copy.deepcopy(p_all)
        q_coord: npt.NDArray[np.float64] = copy.deepcopy(q_all)

    else:
        assert p_view is not None
        assert q_view is not None
        p_coord = copy.deepcopy(p_all[p_view])
        q_coord = copy.deepcopy(q_all[q_view])

    # Recenter to centroid
    p_cent: npt.NDArray[np.float64] = centroid(p_coord)
    q_cent: npt.NDArray[np.float64] = centroid(q_coord)
    p_coord -= p_cent
    q_coord -= q_cent

    # rmsd_method: RmsdCallable
    if (args.add_idx is None) and (args.remove_idx is None) and (args.ignore_Hydrogen is False):
        idx_atom1 = np.arange(
            len(p_all_atoms), dtype=np.int64)
        idx_atom1 = idx_atom1+1

    coord_square: dict[int, float]
    result_rmsd: float
    coord_square, result_rmsd = kabsch_rmsd(
        p_coord, q_coord, list(idx_atom1))

    if __name__ == "__main__":
        print(f"{" RMSD":>5s}", end=" ")
        print(f"{result_rmsd:>10.5f}")

    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(xyz_tmp)

    if len(idx_atom1) == 0:
        raise ValueError(str("The value of idx_atom1 is error"))
    elif len(coord_square) == 0:
        raise ValueError(str("The value of coord_square is error"))
    else:
        return coord_square, result_rmsd
