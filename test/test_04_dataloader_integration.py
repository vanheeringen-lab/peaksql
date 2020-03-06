import unittest
from torch.utils.data import DataLoader

import peaksql
from test.test_02_database import DATABASE_BED


class TestMLIntegration(unittest.TestCase):
    """ A test class to test integration with common machine learning practices """

    def test_401_iterable(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, nr_rand_pos=100, seq_length=3)
        for seq, label in dataset:
            assert seq.shape == (3, 4,)
            assert label.shape == (1,)

    def test_402_Integration_PyTorch_DataLoader(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, nr_rand_pos=100, seq_length=3)
        dataloader = DataLoader(dataset, batch_size=10)
        for seq, label in dataloader:
            assert tuple(seq.shape) == (10, 3, 4,)
            assert tuple(label.shape) == (10, 1,)
