# Should database be able to generate these? db.get_NarrowPeakDataSet(query)
# answer: probably yes, TODO
import numpy as np
from multiprocessing import Lock

from .util import sequence_to_onehot


class NarrowPeakDataSet:
    """
    Dataset ...
    """
    def __init__(self, database, query, size, stride):
        """

        """
        self.size = size
        self.stride = stride
        self.database = database

        database.cursor.execute(query)
        self.fetchall = database.cursor.fetchall()

        # TODO proper test
        assert len(self.fetchall[0]) == 4, "TODO"

        # add desc
        self.peaks = dict()
        for assembly, chrom, chromstart, chromstop in self.fetchall:
            if assembly not in self.peaks:
                self.peaks[assembly] = dict()
            if chrom not in self.peaks[assembly]:
                self.peaks[assembly][chrom] = np.zeros(len(database.fastas[assembly][chrom]), dtype=bool)
            self.peaks[assembly][chrom][chromstart:chromstop] = True

        # add desc
        self.lens, self.order = [0], [None]
        self.locks = {}
        for assembly in self.peaks:
            self.locks[assembly] = Lock()
            for chrom in self.peaks[assembly]:
                self.lens.append((len(self.peaks[assembly][chrom]) - size) // stride)
                self.order.append((assembly, chrom))
        self.lens = np.cumsum(self.lens)

    def __len__(self):
        """

        """
        return self.lens[-1]

    def __getitem__(self, index):
        """

        """
        for i in range(1, len(self.lens)):
            if self.lens[i - 1] <= index < self.lens[i]:
                pos = (index - self.lens[i - 1]) * self.stride
                assembly, chrom = self.order[i]
                chromstart, chromstop = pos, pos + self.size
                with self.locks[assembly]:
                    seq = self.database.fastas[assembly][chrom][chromstart:chromstop]
                return sequence_to_onehot(seq)
        assert False
