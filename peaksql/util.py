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
    return _sequence_to_onehot(np.fromiter(str(sequence).upper(), (np.unicode, 1)))
