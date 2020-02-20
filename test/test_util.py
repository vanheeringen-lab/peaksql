import unittest
import numpy as np

import peaksql.util


SEQ_LEN = 200


class TestBackground(unittest.TestCase):
    """ A test class to test the util of peaksql """

    def test_single_nucleotide_IUPAC_ACGT(self):
        """
        Test whether single-nucleotide IUPAC codes are converted correctly into one-hot encoding.
        """
        sequence = "ACGTACGT"
        onehot = peaksql.util.sequence_to_onehot(sequence)
        true = np.array([[True,  False, False, False],
                         [False, True,  False, False],
                         [False, False, True,  False],
                         [False, False, False, True],
                         [True,  False, False, False],
                         [False, True,  False, False],
                         [False, False, True,  False],
                         [False, False, False, True]], dtype=bool)
        np.testing.assert_array_equal(onehot, true)

    def test_multi_nucleotide_IUPAC_R(self):
        sequence = "".join(["R"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert np.any(onehot[:, 0])
        assert not np.any(onehot[:, 1])
        assert np.any(onehot[:, 2])
        assert not np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_Y(self):
        sequence = "".join(["Y"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert not np.any(onehot[:, 0])
        assert np.any(onehot[:, 1])
        assert not np.any(onehot[:, 2])
        assert np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_S(self):
        sequence = "".join(["S"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert not np.any(onehot[:, 0])
        assert np.any(onehot[:, 1])
        assert np.any(onehot[:, 2])
        assert not np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_W(self):
        sequence = "".join(["W"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert np.any(onehot[:, 0])
        assert not np.any(onehot[:, 1])
        assert not np.any(onehot[:, 2])
        assert np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_K(self):
        sequence = "".join(["K"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert not np.any(onehot[:, 0])
        assert not np.any(onehot[:, 1])
        assert np.any(onehot[:, 2])
        assert np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_M(self):
        sequence = "".join(["M"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert np.any(onehot[:, 0])
        assert np.any(onehot[:, 1])
        assert not np.any(onehot[:, 2])
        assert not np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_B(self):
        sequence = "".join(["B"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert not np.any(onehot[:, 0])
        assert np.any(onehot[:, 1])
        assert np.any(onehot[:, 2])
        assert np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_D(self):
        sequence = "".join(["D"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert np.any(onehot[:, 0])
        assert not np.any(onehot[:, 1])
        assert np.any(onehot[:, 2])
        assert np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_H(self):
        sequence = "".join(["H"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert np.any(onehot[:, 0])
        assert np.any(onehot[:, 1])
        assert not np.any(onehot[:, 2])
        assert np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_V(self):
        sequence = "".join(["V"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert np.any(onehot[:, 0])
        assert np.any(onehot[:, 1])
        assert np.any(onehot[:, 2])
        assert not np.any(onehot[:, 3])

    def test_multi_nucleotide_IUPAC_N(self):
        sequence = "".join(["N"] * SEQ_LEN)
        onehot = peaksql.util.sequence_to_onehot(sequence)
        assert np.any(onehot[:, 0])
        assert np.any(onehot[:, 1])
        assert np.any(onehot[:, 2])
        assert np.any(onehot[:, 3])
