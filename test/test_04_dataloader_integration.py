import unittest
from torch.utils.data import DataLoader
from keras.models import Sequential
from keras.layers import Conv1D, Flatten, Dense

import peaksql
from test.test_02_database import DATABASE_BED


class TestMLIntegration(unittest.TestCase):
    """ A test class to test the peaksql.datasets._BedDataSet family """

    def test_401_iterable(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, nr_rand_pos=100, seq_length=3)
        for seq, label in dataset:
            pass

    def test_402_Integration_PyTorch_DataLoader(self):
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, nr_rand_pos=100, seq_length=3)
        dataloader = DataLoader(dataset, batch_size=8)
        for seq, label in dataloader:
            pass

    def test_403_Integration_Keras_SimpleCNN(self):
        model = Sequential()
        model.add(Conv1D(2, 4, input_shape=(4, 6)))
        model.add(Flatten())
        model.add(Dense(4, activation="sigmoid"))

        model.compile(
            loss="binary_crossentropy", optimizer="adam", metrics=["accuracy"]
        )

        # feed our dataloader to the model and learn
        dataset = peaksql.BedRegionDataSet(DATABASE_BED, nr_rand_pos=100, seq_length=6)
        for seq, label in iter(dataset):
            print(seq, label)
        model.fit_generator(iter(dataset), steps_per_epoch=100)
