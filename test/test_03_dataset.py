import unittest
import numpy as np

import peaksql
from test.test_02_database import DATABASE_BED, DATABASE_NWP


class TestDataBase(unittest.TestCase):
    """ A test class to test the peaksql.datasets """

    positions = np.array(
        [
            [False, True, False, False, False],
            [False, True, True, False, False],
            [False, False, False, False, False],
            [True, True, True, True, True],
        ]
    )

    def test_301_BedDataSet_stride_length(self):
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert len(dataset) == 16

    def test_302_BedDataSet_stride_sequences(self):
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, stride=10)

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
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert all(dataset.any(self.positions) == [True, True, False, True])

    def test_304_Bed_label_all(self):
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert all(dataset.all(self.positions) == [False, False, False, True])

    def test_305_Bed_label_fraction(self):
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, stride=10)
        dataset.ratio = 0.4
        assert all(dataset.fraction(self.positions) == [False, True, False, True])

    def test_306_BedDataSet_array_from_query(self):
        chromstart = 10
        chromend = 20
        query = [(0, 15, 25), (0, 5, 13)]
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert np.all(
            dataset.array_from_query(query, chromstart, chromend)
            == np.array(
                [[True, True, True, False, False, True, True, True, True, True]]
            )
        )

    def test_307_NarrowPeak_array_from_query(self):
        chromstart = 10
        chromend = 20
        query = [(0, 15, 25), (0, 5, 13)]
        dataset = peaksql.NarrowPeakDataSet(DATABASE_NWP, seq_length=10, stride=10)
        assert np.all(
            dataset.array_from_query(query, chromstart, chromend)
            == np.array(
                [[False, False, False, False, False, False, False, False, True, False]]
            )
        )

    def test_308_BedDataSet_random_pos_length(self):
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, nr_rand_pos=20)
        assert len(dataset) == 20

    def test_309_BedDataSet_random_pos_sequences(self):
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, nr_rand_pos=20)
        all_dna = [
            "AAAACCCCGGGGTTTTAAACCCGGGTTTAACCGGTTACGT",
            "TTTTGGGGCCCCAAAATTTGGGCCCAAATTGGCCAATGCA",
            "ATGCGTAGCTGATCGATGCTAGCTAGCTAGCTAGCTAAAA",
            "ATGGTGAATGTGAGTAGTGATGATGAGTGTAGTGAGGGGG",
        ]
        dna_onehot = [peaksql.util.sequence_to_onehot(dna) for dna in all_dna]
        for i, (seq, label) in enumerate(dataset):
            found = False
            for chromosome in range(4):
                for idx in range(0, 30):
                    if np.all(seq == dna_onehot[chromosome][idx : idx + 10]):
                        found = True
            assert found

    def test_310_BedDataSet_random_pos_distribution(self):
        dataset = peaksql.BedDataSet(DATABASE_BED, seq_length=10, nr_rand_pos=100_000)

        # chromosomes are of equal size, so we expect equal nr of positions for each
        un_cumsum = dataset.cumsum - np.roll(dataset.cumsum, shift=1)
        for count in un_cumsum[1:]:
            assert 0.245 <= count / 100_000 <= 0.255
