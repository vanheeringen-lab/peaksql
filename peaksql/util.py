# TODO: random seed + tests
import numba
import numpy as np


@numba.jit(nopython=True, cache=True)
def nuc_to_onehot(nuc):
    """
    Convert a nucleotide to a one hot encoded boolean array, where the indexes respectively
    correspond to A, C, G, T.

    Accepts all IUPAC nucleotide codes, and picks a random option from the possible nucleotides.
    """
    # first try single-nucleotide codes
    if nuc == 'A':
        return np.array([1, 0, 0, 0], dtype=numba.boolean)
    if nuc == 'C':
        return np.array([0, 1, 0, 0], dtype=numba.boolean)
    if nuc == 'G':
        return np.array([0, 0, 1, 0], dtype=numba.boolean)
    if nuc == 'T':
        return np.array([0, 0, 0, 1], dtype=numba.boolean)

    # if that doesn't work try multiple-nucleotide codes
    # first set potential indexes
    if nuc == 'N':
        idx = np.array([0, 1, 2, 3])
    elif nuc == 'R':
        idx = np.array([0, 2])
    elif nuc == 'Y':
        idx = np.array([1, 3])
    elif nuc == 'S':
        idx = np.array([1, 2])
    elif nuc == 'W':
        idx = np.array([0, 3])
    elif nuc == 'K':
        idx = np.array([2, 3])
    elif nuc == 'M':
        idx = np.array([0, 1])
    elif nuc == 'B':
        idx = np.array([1, 2, 3])
    elif nuc == 'D':
        idx = np.array([0, 2, 3])
    elif nuc == 'H':
        idx = np.array([0, 1, 3])
    elif nuc == 'V':
        idx = np.array([0, 1, 2])
    else:
        raise ValueError("Only IUPAC nucleotide codes are accepted.")

    # now make an empty encoding, and set a random index to True
    onehot = np.array([0, 0, 0, 0], dtype=numba.boolean)
    onehot[np.random.choice(idx)] = 1

    return onehot


@numba.jit(nopython=True, cache=True)
def _sequence_to_onehot(sequence):
    onehot = np.zeros((len(sequence), 4), dtype=numba.boolean)
    for i, nuc in enumerate(sequence):
        onehot[i] = nuc_to_onehot(nuc)

    return onehot

@profile
def sequence_to_onehot(sequence):
    """
    Convert a sequence of length n to one-hot encoding of shape (4 x n).
    """
    return _sequence_to_onehot(np.fromiter(str(sequence).upper(), (np.unicode, 1))).T


@numba.jit(nopython=True, cache=True)
def binary_search(index, lens):
    """
    Does a binary search to not find the value in a list, but the index where the value is in
    between two values (returns the higher index of the two).

    index: 5
    lens: [0, 4, 6, 8, 22]
    returns: 2
    """
    left, right = 0, len(lens)

    while left <= right:
        mid = (left + right) // 2

        if lens[mid - 1] > index:
            right = mid

        elif lens[mid] <= index:
            left = mid

        else:
            assert lens[mid - 1] <= index < lens[mid]
            return mid

    assert False
