import unittest
import numpy as np

import peaksql.util


SEQ_LEN = 200


class TestUtil(unittest.TestCase):
    """ A test class to test the util of peaksql """

    def test_101_IUPAC_A_to_idx(self):
        onehot_idx = peaksql.util._nuc_to_onehot_idx.py_func(ord("A"))
        assert onehot_idx == 0

    def test_102_IUPAC_C_to_idx(self):
        onehot_idx = peaksql.util._nuc_to_onehot_idx.py_func(ord("C"))
        assert onehot_idx == 1

    def test_103_IUPAC_G_to_idx(self):
        onehot_idx = peaksql.util._nuc_to_onehot_idx.py_func(ord("G"))
        assert onehot_idx == 2

    def test_104_IUPAC_T_to_idx(self):
        onehot_idx = peaksql.util._nuc_to_onehot_idx.py_func(ord("T"))
        assert onehot_idx == 3

    def test_105_IUPAC_R_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("R")) for _ in range(SEQ_LEN)]
        )
        assert np.any(onehot_idx == 0)
        assert not np.any(onehot_idx == 1)
        assert np.any(onehot_idx == 2)
        assert not np.any(onehot_idx == 3)

    def test_106_IUPAC_Y_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("Y")) for _ in range(SEQ_LEN)]
        )
        assert not np.any(onehot_idx == 0)
        assert np.any(onehot_idx == 1)
        assert not np.any(onehot_idx == 2)
        assert np.any(onehot_idx == 3)

    def test_107_IUPAC_S_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("S")) for _ in range(SEQ_LEN)]
        )
        assert not np.any(onehot_idx == 0)
        assert np.any(onehot_idx == 1)
        assert np.any(onehot_idx == 2)
        assert not np.any(onehot_idx == 3)

    def test_108_IUPAC_W_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("W")) for _ in range(SEQ_LEN)]
        )
        assert np.any(onehot_idx == 0)
        assert not np.any(onehot_idx == 1)
        assert not np.any(onehot_idx == 2)
        assert np.any(onehot_idx == 3)

    def test_109_IUPAC_K_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("K")) for _ in range(SEQ_LEN)]
        )
        assert not np.any(onehot_idx == 0)
        assert not np.any(onehot_idx == 1)
        assert np.any(onehot_idx == 2)
        assert np.any(onehot_idx == 3)

    def test_110_IUPAC_M_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("M")) for _ in range(SEQ_LEN)]
        )
        assert np.any(onehot_idx == 0)
        assert np.any(onehot_idx == 1)
        assert not np.any(onehot_idx == 2)
        assert not np.any(onehot_idx == 3)

    def test_111_IUPAC_B_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("B")) for _ in range(SEQ_LEN)]
        )
        assert not np.any(onehot_idx == 0)
        assert np.any(onehot_idx == 1)
        assert np.any(onehot_idx == 2)
        assert np.any(onehot_idx == 3)

    def test_112_IUPAC_D_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("D")) for _ in range(SEQ_LEN)]
        )
        assert np.any(onehot_idx == 0)
        assert not np.any(onehot_idx == 1)
        assert np.any(onehot_idx == 2)
        assert np.any(onehot_idx == 3)

    def test_113_IUPAC_H_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("H")) for _ in range(SEQ_LEN)]
        )
        assert np.any(onehot_idx == 0)
        assert np.any(onehot_idx == 1)
        assert not np.any(onehot_idx == 2)
        assert np.any(onehot_idx == 3)

    def test_114_IUPAC_V_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("V")) for _ in range(SEQ_LEN)]
        )
        assert np.any(onehot_idx == 0)
        assert np.any(onehot_idx == 1)
        assert np.any(onehot_idx == 2)
        assert not np.any(onehot_idx == 3)

    def test_115_IUPAC_N_to_idx(self):
        onehot_idx = np.array(
            [peaksql.util._nuc_to_onehot_idx.py_func(ord("N")) for _ in range(SEQ_LEN)]
        )
        assert np.any(onehot_idx == 0)
        assert np.any(onehot_idx == 1)
        assert np.any(onehot_idx == 2)
        assert np.any(onehot_idx == 3)

    def test_116_IUPAC_wrong(self):
        self.assertRaises(ValueError, peaksql.util._nuc_to_onehot_idx.py_func, ord("Q"))

    def test_117_sequence_to_onehot(self):
        """
        Test whether single-nucleotide IUPAC codes are converted correctly into one-hot
        encoding.
        """
        sequence = "ACGTACGT"
        onehot = peaksql.util.sequence_to_onehot(sequence)
        true = np.array(
            [
                [True, False, False, False],
                [False, True, False, False],
                [False, False, True, False],
                [False, False, False, True],
                [True, False, False, False],
                [False, True, False, False],
                [False, False, True, False],
                [False, False, False, True],
            ],
            dtype=bool,
        )
        np.testing.assert_array_equal(onehot, true)

    def test_118_binary_search_out_of_bound(self):
        haystack = np.array([0, 5, 10, 12, 22])

        # call each function twice, once jit compiled and once not compiled for coverage
        self.assertRaises(ValueError, peaksql.util.binary_search, -1, haystack)
        self.assertRaises(ValueError, peaksql.util.binary_search.py_func, -1, haystack)

        self.assertRaises(ValueError, peaksql.util.binary_search, 22, haystack)
        self.assertRaises(ValueError, peaksql.util.binary_search.py_func, 22, haystack)

    def test_119_binary_search(self):
        haystack = np.array([0, 5, 10, 12, 22])

        # call each function twice, once jit compiled and once not compiled for coverage
        assert peaksql.util.binary_search(0, haystack) == 1
        assert peaksql.util.binary_search.py_func(0, haystack) == 1

        assert peaksql.util.binary_search(1, haystack) == 1
        assert peaksql.util.binary_search.py_func(1, haystack) == 1

        assert peaksql.util.binary_search(14, haystack) == 4
        assert peaksql.util.binary_search.py_func(14, haystack) == 4
