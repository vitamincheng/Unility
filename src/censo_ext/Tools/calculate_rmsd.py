#!/usr/bin/env python
import argparse
__doc__ = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation.

For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

# from icecream import ic
import copy
import numpy as np
import numpy.typing as npt
from pathlib import Path
import censo_ext.Tools.Parameter as Parameter
from censo_ext.Tools.xyzfile import GeometryXYZs


NAMES_ELEMENT = {value: key for key, value in Parameter.ELEMENT_NAMES.items()}


def atom2str(atom: int) -> str:
    """
    Convert atom type from integer to string.

    Args:
        atom (int): Atom type as integer.

    Returns:
        str: Atom type as string.
    """
    return Parameter.ELEMENT_NAMES[atom]


def atom2int(atom: str) -> int:
    """
    Convert atom type from string to integer.

    Args:
        atom (str): Atom type as string.

    Returns:
        int: Atom type as integer.
    """
    atom = atom.capitalize().strip()
    return NAMES_ELEMENT[atom]


def rmsd(P: npt.NDArray[np.float64], Q: npt.NDArray[np.float64], idx_atom: list[int], **kwargs) -> tuple[dict[int, float], float]:
    """
    Calculate Root-mean-square deviation from two sets of vectors.

    Args:
        P(ndarray): (N,D) matrix, where N is points and D is dimension.
        Q(ndarray): (N,D) matrix, where N is points and D is dimension.
        idx_atom(list): List of atom indices to consider in the calculation.
        **kwargs: Additional keyword arguments (not used in current implementation).

    Returns:
        Tuple[dict,float]: A tuple containing:
            - atom_square (dict): Dictionary mapping atom indices to their squared
              coordinate differences.
            - rmsd (float): Root-mean-square deviation between the two vectors.

    Example:
        >>> P = np.array([[1, 2], [3, 4]])
        >>> Q = np.array([[2, 3], [4, 5]])
        >>> idx_atom = [0, 1]
        >>> atom_square, rmsd_value = rmsd(P, Q, idx_atom)
        >>> print(f"Atom square differences: {atom_square}")
        >>> print(f"RMSD: {rmsd_value}")
    """

    diff: npt.NDArray[np.float64] = P - Q
    idx_atomSquare: dict[int, float] = {}
    coord_square_total: float = 0
    for idx0, x in enumerate(idx_atom):
        coord_square: float = float((diff[idx0]*diff[idx0]).sum())
        # ic(coord_square)
        if __name__ == "__main__":
            print(f"{x:>5}", end=" ")
            print(f"{coord_square:>10.5f}")
        idx_atomSquare[x] = coord_square
        coord_square_total += coord_square
    return idx_atomSquare, float(np.sqrt(coord_square_total / P.shape[0]))


def kabsch_rmsd(P: npt.NDArray[np.float64], Q: npt.NDArray[np.float64], idx_atom1: list[int],
                translate: bool = False) -> tuple[dict[int, float], float]:
    """
    Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.

    Args:
        P (npt.NDArray[np.flaot64]): (N,D) matrix, where N is points and D is dimension.
        Q (npt.NDArray[np.flaot64]): (N,D) matrix, where N is points and D is dimension.
        idx_atom1 (list[int]): List of atom indices to consider.
        translate (bool, optional): Use centroids to translate vector P and Q unto each other. Defaults to False.

    Returns:
        tuple[dict[int,float],float]: Tuple of atom-wise squared differences and RMSD value.
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

    Args:
        P (npt.NDArray[np.float64]): (N,D) matrix, where N is points and D is dimension.
        Q (npt.NDArray[np.float64]): (N,D) matrix, where N is points and D is dimension.

    Returns:
        npt.NDArray[np.float64]: Rotated matrix P.
    """
    U: npt.NDArray[np.float64] = kabsch(P, Q)

    # Rotate P
    return np.dot(P, U)


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

    Args:
        P (npt.NDArray[np.float64]): (N,D) matrix, where N is points and D is dimension.
        Q (npt.NDArray[np.float64]): (N,D) matrix, where N is points and D is dimension.

    Returns:
        npt.NDArray[np.float64]: Rotation matrix (D,D)
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
    return np.dot(V, W)


def centroid(X: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """
    Centroid is the mean position of all the points in all of the coordinate
    directions, from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid

    C = sum(X)/len(X)

    Args:
        X (npt.NDArray[np.float64]): (N,D) matrix, where N is points and D is dimension.

    Returns:
        npt.NDArray[np.float64]: Centroid.
    """
    return X.mean(axis=0)


def get_Coordinates(xyzFile, idx0) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    """Read xyz file to extract atom types and coordinates.

    Args:
        xyzFile: GeometryXYZs object containing molecular structures.
        idx0(int): Index of the structure in xyzfile (zero-based indexing).

    Returns:
        tuple[npt.NDArray[npt.NDArray[np.int64]],npt.NDArray[np.float64]]: 
        A tuple containing:
        - element: Array of atom types as integers
        - V: Array of atomic coordinates as floats

    Example:
        >>> element, coordinates = get_Coordinates(xyzfile, 0)
        >>> print(element)
        [6 8 1 1]
        >>> print(coordinates)
        [[ 0.0  0.0  0.0]
         [ 0.0  0.0  1.0]
         [ 1.0  0.0  0.0]
         [ 0.0  1.0  0.0]]
    """

    idx_Names: dict[int, str] = xyzFile.Sts[idx0].names
    element: npt.NDArray[np.int64] = np.array(
        [atom2int(atom) for atom in idx_Names.values()])
    V: npt.NDArray[np.float64] = np.array(
        xyzFile.Sts[idx0].coord, dtype=np.float64)
    return element, V


def cal_RMSD_xyz(xyzFile: GeometryXYZs, idx_p: int, idx_q: int, args: argparse.Namespace) -> tuple[dict[int, float], float]:
    """
    Read xyz file and calculate RMSD between two structures.

    This function reads two structures from an XYZ file, processes them according
    to the provided arguments (such as ignoring hydrogen atoms, removing or adding 
    specific atom indices, breaking bonds), and calculates the RMSD using the Kabsch algorithm.

    Args:
        xyzFile(GeometryXYZs): The XYZ file object containing molecular structures.
        idx_p(int): Index of the first structure (1-based indexing).
        idx_q(int): Index of the second structure (1-based indexing).
        args(argparse.Namespace): Command-line arguments controlling the RMSD calculation.

    Returns:
        tuple[dict[int,float],float]: A tuple containing:
            - dict[int,float]: Dictionary mapping atom indices to their squared differences.
            - float: The final RMSD value.

    Raises:
        ValueError: If structures have different sizes or if idx_atom1 is empty after processing.
    """

    xyz_tmp: Path = Path(".tmp.xyz")
    idx_p -= 1
    idx_q -= 1
    p_all_atoms: npt.NDArray[np.int64]
    q_all_atoms: npt.NDArray[np.int64]
    p_all: npt.NDArray[np.float64]
    q_all: npt.NDArray[np.float64]
    p_all_atoms, p_all = get_Coordinates(xyzFile, idx_p)
    q_all_atoms, q_all = get_Coordinates(xyzFile, idx_q)

    if p_all.shape[0] != q_all.shape[0]:
        raise ValueError("error: Structures not same size")

    # Initialize atom indices
    idx_atom1: npt.NDArray[np.int64] = np.array([], dtype=np.int64)
    index: set[int] | list[int] | npt.NDArray[np.int64]
    p_view: None | npt.NDArray[np.int64] = None
    q_view: None | npt.NDArray[np.int64] = None

    # Handle hydrogen removal
    if args.ignore_Hydrogen:
        assert type(p_all_atoms[0]) is not str
        assert type(q_all_atoms[0]) is not str
        for idx in range(len(p_all_atoms)):
            if p_all_atoms[idx] != 1:
                idx_atom1 = np.append(idx_atom1, [idx+1])

        p_view = np.where(p_all_atoms != 1)  # type: ignore
        q_view = np.where(q_all_atoms != 1)  # type: ignore

    # Handle bond breaking
    if args.bond_broken:
        if args.ignore_Hydrogen:
            xyzFile.set_filename(xyz_tmp)
            xyzFile.method_save_xyz([idx_p])
            args_x: dict = {"file": xyz_tmp, "bond_broken": [*args.bond_broken],
                            "print": False, "debug": False}
            from censo_ext.Tools.topo import Topo
            Sts_topo: Topo = Topo(args_x["file"])
            idx_atom1: npt.NDArray[np.int64] = np.array(Sts_topo.method_broken_bond(
                argparse.Namespace(**args_x)))
            index = idx_atom1-1
            p_view, q_view = index, index

        else:
            raise ValueError(" Only support under ignore Hydrogen condition ")
    else:
        pass

    # Handle index removal
    if args.remove_idx:
        if not args.ignore_Hydrogen:
            idx_atom1 = np.arange(1, len(p_all_atoms)+1)

        idx_atom1 = np.setdiff1d(idx_atom1, args.remove_idx)

        args.remove_idx = np.array(args.remove_idx)-1
        index = idx_atom1-1

        p_view, q_view = index, index

    # Handle index addition
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

    # Final index handling
    if (args.add_idx is None) and (args.remove_idx is None) and (not args.ignore_Hydrogen):
        idx_atom1 = np.arange(1, len(p_all_atoms)+1)

    idx_coordSquare: dict[int, float]
    res_rmsd: float
    idx_coordSquare, res_rmsd = kabsch_rmsd(p_coord, q_coord, list(idx_atom1))

    if __name__ == "__main__":
        print(f"{" RMSD":>5s}", end=" ")
        print(f"{res_rmsd:>10.5f}")

    from censo_ext.Tools.utility import delete_all_files
    delete_all_files(xyz_tmp)

    if len(idx_atom1) == 0:
        raise ValueError("The value of idx_atom1 is error")
    elif len(idx_coordSquare) == 0:
        raise ValueError("The value of coord_square is error")
    else:
        return idx_coordSquare, res_rmsd
