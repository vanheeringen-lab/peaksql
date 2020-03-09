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
        query = [[0, 0, 15, 25], [1, 0, 5, 13]]
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, stride=10)
        assert np.all(
            dataset._array_from_query(query, chromid, chromstart, chromend)
            == np.array(
                [[True, True, True, False, False, False, False, False, False, False]]
            )
        )

    def test_310_NarrowPeak_array_from_query(self):
        chromstart = 10
        chromend = 20
        chromid = 1
        query = [[0, 0, 15, 25], [1, 0, 5, 13]]
        dataset = peaksql.NarrowPeakDataSet(DATABASE_NWP, seq_length=10, stride=10)
        assert np.all(
            dataset._array_from_query(query, chromid, chromstart, chromend)
            == np.array(
                [[False, False, False, False, False, False, False, False, True, False]]
            )
        )

    def test_311_BedRegionDataSet_random_pos_length(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, nr_rand_pos=20)
        assert len(dataset) == 20

    def test_312_BedRegionDataSet_random_pos_sequences(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, seq_length=10, nr_rand_pos=20)
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

    def test_313_BedRegionDataSet_random_pos_distribution(self):
        dataset = peaksql.BedRegionDataSet(
            DATABASE_BED, seq_length=10, nr_rand_pos=100_000
        )

        # chromosomes are of equal size, so we expect equal nr of positions for each
        un_cumsum = dataset.cumsum - np.roll(dataset.cumsum, shift=1)
        for count in un_cumsum[1:]:
            assert 0.245 <= count / 100_000 <= 0.255

    def test_314_BedDataSet_label_func_parsing(self):
        with self.assertRaises(ValueError):
            peaksql.BedRegionDataSet(
                DATABASE_BED, seq_length=10, nr_rand_pos=10, label_func="inner_any"
            )
        with self.assertRaises(ValueError):
            peaksql.BedRegionDataSet(
                DATABASE_BED, seq_length=10, nr_rand_pos=10, label_func="inner_all"
            )
        with self.assertRaises(ValueError):
            peaksql.BedRegionDataSet(
                DATABASE_BED,
                seq_length=10,
                nr_rand_pos=10,
                label_func="inner_fraction",
                fraction=0.8,
            )
        with self.assertRaises(ValueError):
            peaksql.BedRegionDataSet(
                DATABASE_BED, seq_length=10, nr_rand_pos=10, label_func="fraction"
            )
