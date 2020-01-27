import numba
import numpy as np


@numba.jit(nopython=True)
def nuc_to_onehot(nuc):
    if nuc == 'A':
        return [1.00, 0.00, 0.00, 0.00]
    if nuc == 'C':
        return [0.00, 1.00, 0.00, 0.00]
    if nuc == 'G':
        return [0.00, 0.00, 1.00, 0.00]
    if nuc == 'T':
        return [0.00, 0.00, 0.00, 1.00]
    if nuc == 'N':
        return [0.25, 0.25, 0.25, 0.25]
    raise ValueError


@numba.jit(nopython=True)
def _sequence_to_onehot(sequence):
    onehot = np.zeros((len(sequence), 4))
    for i, nuc in enumerate(sequence):
        onehot[i] = nuc_to_onehot(sequence[i])

    return onehot


def sequence_to_onehot(sequence):
    """
    Convert a sequence of length n to one-hot encoding of shape (4 x n).
    """
    return _sequence_to_onehot(np.fromiter(str(sequence).upper(), (np.unicode, 1))).T


@numba.jit(nopython=True)
def _binary_search(index, lens):
    """
    Does a binary search to not find the value in a list, but the index where the value is in
    between two values.
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


@numba.jit(nopython=True)
def jit_any(arr):
    """
    Jit-compiled numpy.any function.
    """
    return np.any(arr)
