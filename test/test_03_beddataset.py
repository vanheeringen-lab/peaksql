import unittest
import numpy as np

import peaksql
from test.test_02_database import DATABASE_BED, DATABASE_NWP


class TestDataBase(unittest.TestCase):
    """ A test class to test the peaksql.datasets._BedDataSet family """

    positions = np.array(
        [
            [False, True, False, False, False],
            [False, True, True, False, False],
            [False, False, False, False, False],
            [True, True, True, True, True],
        ]
    )

    def test_301_BedRegionDataSet_stride_length(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert len(dataset) == 16

    def test_302_BedRegionDataSet_stride_sequences(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)

        all_dna = (
            "AAAACCCCGGGGTTTTAAACCCGGGTTTAACCGGTTACGT"
            + "TTTTGGGGCCCCAAAATTTGGGCCCAAATTGGCCAATGCA"
            + "ATGCGTAGCTGATCGATGCTAGCTAGCTAGCTAGCTAAAA"
            + "ATGGTGAATGTGAGTAGTGATGATGAGTGTAGTGAGGGGG"
        )

        dna_strided = [all_dna[i : i + 10] for i in range(0, len(all_dna), 10)]
        dna_onehot = [peaksql.util.sequence_to_onehot(dna) for dna in dna_strided]
        for seq, label in dataset:
            assert np.sum(np.all(seq == potential_seq) for potential_seq in dna_onehot)

    def test_303_Bed_label_any(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert all(dataset.label_any(self.positions) == [True, True, False, True])

    def test_304_Bed_label_all(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert all(dataset.label_all(self.positions) == [False, False, False, True])

    def test_305_Bed_label_fraction(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        dataset.fraction = 0.4
        assert all(dataset.label_fraction(self.positions) == [False, True, False, True])

    def test_306_Bed_label_inner_any(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        dataset.inner_range = 1
        assert all(dataset.label_any(self.positions) == [True, True, False, True])

    def test_307_Bed_label_inner_all(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        dataset.inner_range = 1
        assert all(
            dataset.label_inner_all(self.positions) == [False, False, False, True]
        )

    def test_308_Bed_label_inner_fraction(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        dataset.inner_range = 1
        dataset.fraction = 2 / 3
        assert all(
            dataset.label_inner_fraction(self.positions) == [False, True, False, True]
        )

    def test_309_BedRegionDataSet_array_from_query(self):
        chromstart = 10
        chromend = 20
        chromid = 1
        query = [[0, None, 15, 25], [1, None, 5, 13]]
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert np.all(
            dataset.array_from_query(query, chromid, chromstart, chromend)
            == np.array(
                [[True, True, True, False, False, False, False, False, False, False]]
            )
        )

    def test_309_NarrowPeak_array_from_query(self):
        chromstart = 10
        chromend = 20
        chromid = 1
        query = [[0, None, 15, 25], [1, None, 5, 13]]
        dataset = peaksql.NarrowPeakDataSet(DATABASE_NWP, seq_length=10, stride=10)
        assert np.all(
            dataset.array_from_query(query, chromid, chromstart, chromend)
            == np.array(
                [[False, False, False, False, False, False, False, False, True, False]]
            )
        )
