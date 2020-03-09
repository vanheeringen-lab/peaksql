# TODO: random seed + tests
import numba
import numpy as np


@numba.jit(nopython=True, cache=True)
def _nuc_to_onehot_idx(nuc: int) -> int:
    """
    Convert a nucleotide to a one hot index, where the indexes 0, 1, 2, 3 respectively
    correspond to A, C, G, T.

    Accepts all IUPAC nucleotide codes, and picks a random option from the possible
    nucleotides.
    """
    # first try single-nucleotide codes
    # start with A & T, since genomes generally have higher AT content than GC content
    if nuc == 65:  # A
        return 0
    if nuc == 84:  # T
        return 3
    if nuc == 67:  # C
        return 1
    if nuc == 71:  # G
        return 2

    # if that doesn't work try multiple-nucleotide codes
    # first set potential indexes
    if nuc == 78:  # N
        idx = np.array([0, 1, 2, 3])
    elif nuc == 82:  # R
        idx = np.array([0, 2])
    elif nuc == 89:  # Y
        idx = np.array([1, 3])
    elif nuc == 83:  # S
        idx = np.array([1, 2])
    elif nuc == 87:  # W
        idx = np.array([0, 3])
    elif nuc == 75:  # K
        idx = np.array([2, 3])
    elif nuc == 77:  # M
        idx = np.array([0, 1])
    elif nuc == 66:  # B
        idx = np.array([1, 2, 3])
    elif nuc == 68:  # D
        idx = np.array([0, 2, 3])
    elif nuc == 72:  # H
        idx = np.array([0, 1, 3])
    elif nuc == 86:  # V
        idx = np.array([0, 1, 2])
    else:
        raise ValueError("Only IUPAC nucleotide codes are accepted.")

    # now return any of the possible indexes
    return np.random.choice(idx)


@numba.jit(nopython=True, cache=True)
def _sequence_to_onehot(sequence: np.ndarray) -> np.ndarray:
    onehot = np.zeros((len(sequence), 4), dtype=numba.boolean)
    for i, nuc in enumerate(sequence):
        onehot[i, _nuc_to_onehot_idx(nuc)] = True

    return onehot


def sequence_to_onehot(sequence) -> np.ndarray:
    """
    Convert a sequence of length n to a one-hot encoded array of shape (n x 4).

    The nucleotides A, C, G, T respectively correspond to indices 0, 1, 2, 3.
    """
    return _sequence_to_onehot(str(sequence).upper().encode("utf-8"))


@numba.jit(nopython=True, cache=True)
def binary_search(index: int, lens: np.ndarray) -> int:
    """
    Does a binary search to not find the value in a list, but the index where the value
    is in between two values (returns the higher index of the two).

    index: 5
    lens: [0, 4, 6, 8, 22]
    returns: 2
    """
    if index < lens[0] or index >= lens[-1]:
        raise ValueError("Invalid index")

    # FIXME: this can hang in infinite loop, and because of numba becomes unresponsive.
    left, right = 0, len(lens)

    while left <= right:
        mid = (left + right) // 2

        if lens[mid - 1] > index:
            right = mid

        elif lens[mid] <= index:
            left = mid

        else:
            return mid

    assert False
