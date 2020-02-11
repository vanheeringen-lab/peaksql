import numpy as np
import multiprocessing
from abc import ABC, abstractmethod

import pyfaidx
import sqlite3

from database import DataBase
from util import sequence_to_onehot, _binary_search, jit_any, get_label_count


class DataSet(ABC):
    """

    """
    def __init__(self, database: str, query: str = "", seq_length: int = 200, **kwargs: int):
        if ('stride' in kwargs) == ('nr_rand_pos' in kwargs):  # xor
            raise ValueError("choose a stride or a number of random positions")

        #
        self.database = DataBase(database)

        # sql lookup
        query = self.SELECT + query
        self.database.cursor.execute(query)
        self.fetchall = self.database.cursor.fetchall()


        #
        self.seq_length = seq_length
        self.fastas = {}
        if "stride" in kwargs:
            self.stride = kwargs['stride']
            self.chromosomes, self.cumsum, self.positions = \
                self.get_strided_positions(self.seq_length, self.stride)
        if "nr_rand_pos" in kwargs:
            self.nr_rand_pos = kwargs['nr_rand_pos']
            self.chromosomes, self.cumsum, self.positions = \
                self.get_random_positions(self.seq_length, self.nr_rand_pos)

    def __len__(self) -> int:
        """
        Return the number of indices this dataset contains.
        """
        return self.cumsum[-1]

    @abstractmethod
    def __getitem__(self, index: int):
        pass

    def _get_process(self):
        """
        PyFaidx is not multiprocessing safe when reading from fasta index. However if we start a
        Fasta class for each process, the file pointers don't get mangled, and we can use
        multiprocessing. This allows us to "stream" our data parallel in e.g. a Pytorch dataloader.
        """
        process = multiprocessing.current_process().name
        if process not in self.fastas:
            # start a new connection
            conn = sqlite3.connect(self.database.db)
            cursor = conn.cursor()

            cursor.execute("SELECT Assembly, AbsPath FROM Assembly")
            self.fastas[process] = {
                assembly: pyfaidx.Fasta(abspath)
                for assembly, abspath in cursor.fetchall()
            }
        return process

    def get_strided_positions(self, seq_length: int, stride: int):
        """
        Calculate a map that connects __getitem__ indices to (assembly, chrom, chromstart) triplet.
        The positions are sampled accross the query with an even stride.

        The first return value is a list of (assembly, chrom) pairs, the second list consists of
        the cumulative sum of the number of indices that belong to this assembly, chrom pair, and
        the third list contains all chromstarts of the sequences. This allows for a decently fast
        and memory-efficient lookup of genomic positions corresponding to an index.
        """
        combis = [(None, None)] + list({(assembly, chrom) for assembly, chrom, *_ in self.fetchall})

        counts = [0]
        startpos = [np.array([])]
        for assembly, chrom in combis[1:]:
            counts.append(len(self.database.fastas[assembly][chrom]) // stride)
            startpos.append(np.arange(0, len(self.database.fastas[assembly][chrom]) - seq_length, seq_length))

        cumsum = np.cumsum(counts)

        return combis, cumsum, startpos

    def get_random_positions(self, seq_length, nr_rand_pos):
        """
        Calculate a map that connects __getitem__ indices to (assembly, chrom, chromstart) triplet.
        The positions are sampled accross the query randomly, but proportional to the size of each
        chromosome.

        The first return value is a list of (assembly, chrom) pairs, the second list consists of
        the cumulative sum of the number of indices that belong to this assembly, chrom pair, and
        the third list contains all chromstarts of the sequences. This allows for a decently fast
        and memory-efficient lookup of genomic positions corresponding to an index.
        """
        combis = [(None, None)] + list({(assembly, chrom) for assembly, chrom, *_ in self.fetchall})

        sizes = []
        for assembly, chrom in combis[1:]:
            sizes.append(len(self.database.fastas[assembly][chrom]))

        distribution = np.random.choice(range(len(sizes)), size=nr_rand_pos, p=np.array(sizes) / np.sum(sizes))

        counts = [0]
        startpos = [np.array([])]
        for i, (assembly, chrom) in enumerate(combis[1:]):
            count = np.sum(distribution == i)
            counts.append(count)
            startpos.append(np.random.randint(0, len(self.database.fastas[assembly][chrom]) - seq_length, size=count))

        cumsum = np.cumsum(counts)

        return combis, cumsum, startpos

    def _index_to_site(self, index: int) -> (str, str, int, int):
        """
        Convert the index of self.__getitem__ to a tuple of (assembly, chrom, chromstart, chromstop)

        Uses binary search for fast retrieval.
        """
        index_bs = _binary_search(index, self.cumsum)

        assembly, chrom = self.chromosomes[index_bs]
        chromstart = self.positions[index_bs][index - self.cumsum[index_bs - 1]]
        chromstop = chromstart + self.seq_length

        return assembly, chrom, chromstart, chromstop


class BedDataSet(DataSet):
    SELECT = "SELECT Assembly, Chromosome, ChromStart, ChromStop FROM Bed " \
             "INNER JOIN Chromosome Chr ON Bed.ChromosomeId = Chr.ChromosomeId " \
             "INNER JOIN Condition Con  ON Bed.ConditionId  = Con.ConditionId " \
             "INNER JOIN Assembly Ass   ON Chr.AssemblyId   = Ass.AssemblyId "

    def __init__(self, database: str, query: str = "", seq_length: int = 200, **kwargs: int):
        DataSet.__init__(self, database, query, seq_length, **kwargs)
        self.peaks = self.get_peak_locations()

    def __getitem__(self, index: int) -> (np.array, int):
        """
        Return the sequence in one-hot encoding and the label of the corresponding index.
        """
        process = self._get_process()
        assembly, chrom, chromstart, chromstop = self._index_to_site(index)

        # get the sequence
        seq = self.fastas[process][assembly][chrom][chromstart:chromstop]
        seq = sequence_to_onehot(seq)

        # get the label
        label = jit_any(self.peaks[assembly][chrom][chromstart:chromstop])

        return seq, label

    def get_peak_locations(self) -> dict:
        """
        Map of all assemblies and genomic coordinates and whether or not it falls within a peak.

        :return
        {
        assembly1:
            {chrom1: np.array([False, False, True, ..., True], dtype=bool),
             chrom2: np.array([False, True, False, ..., True], dtype=bool),
             ...
        }
        """
        peaks = dict()
        for assembly, chrom, chromstart, chromstop in self.fetchall:
            if assembly not in peaks:
                peaks[assembly] = dict()
            if chrom not in peaks[assembly]:
                peaks[assembly][chrom] = np.zeros(
                    len(self.database.fastas[assembly][chrom]),
                    dtype=bool
                )
            peaks[assembly][chrom][chromstart:chromstop] = True

        return peaks

    def get_label_occurence(self) -> (int, int):
        """
        Get the number of occurrences of each label.

        First returns the number of negatives (non-peak), and then the number of positives (peak).
        """
        peaks = 0
        for i, (assembly, chrom) in enumerate(self.chromosomes[1:]):
            for chromstart in self.positions[i + 1]:
                peaks += jit_any(self.peaks[assembly][chrom][chromstart:chromstart+self.seq_length])

        return len(self) - peaks, peaks


class NarrowPeakDataSet(DataSet, BedDataSet):
    SELECT = "SELECT Assembly, Chromosome, ChromStart, Peak FROM Bed " \
             "INNER JOIN Chromosome Chr ON Bed.ChromosomeId = Chr.ChromosomeId " \
             "INNER JOIN Condition Con  ON Bed.ConditionId  = Con.ConditionId " \
             "INNER JOIN Assembly Ass   ON Chr.AssemblyId   = Ass.AssemblyId "



class BamDataSet(DataSet):
    def __init__(self):
        raise NotImplementedError


class WigDataSet(DataSet):
    def __init__(self):
        raise NotImplementedError

import time
now = time.time()
dataset = BedDataSet("/vol/PeakSQL/peakSQL.sqlite", stride=200)
print(dataset.get_label_occurence())
print(time.time() - now)
