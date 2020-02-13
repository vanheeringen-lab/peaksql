import numpy as np
import multiprocessing
from abc import ABC, abstractmethod

import pyfaidx
import sqlite3

from database import DataBase
import util


class DataSet(ABC):
    """

    """
    def __init__(self, database: str, where: str = "", seq_length: int = 200, **kwargs: int):
        if ('stride' in kwargs) == ('nr_rand_pos' in kwargs):  # xor
            raise ValueError("choose a stride or a number of random positions")

        #
        self.database = DataBase(database)

        # sql lookup
        self.WHERE = where
        query = self.SELECT + self.FROM + self.WHERE
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

    def __getitem__(self, index: int) -> (np.array, int):
        """
        Return the sequence in one-hot encoding and the label of the corresponding index.
        """
        process = self._get_process()
        assembly, chrom, chromstart, chromend = self._index_to_site(index)

        # get the sequence, label and condition
        seq = self.get_onehot_sequence(assembly, chrom, chromstart, chromend, process)
        label = self.get_label(assembly, chrom, chromstart, chromend, process)
        condition = self.get_condition(assembly, chrom, chromstart, chromend, process)

        return seq, label, condition

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

    def _index_to_site(self, index: int) -> (str, str, int, int):
        """
        Convert the index of self.__getitem__ to a tuple of (assembly, chrom, chromstart, chromend)

        Uses binary search for fast retrieval.
        """
        index_bs = util.binary_search(index, self.cumsum)

        assembly, chrom = self.chromosomes[index_bs]
        chromstart = self.positions[index_bs][index - self.cumsum[index_bs - 1]]
        chromend = chromstart + self.seq_length

        return assembly, chrom, chromstart, chromend

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

        # distribute the positions over the chromosomes
        sizes = []
        for assembly, chrom in combis[1:]:
            sizes.append(len(self.database.fastas[assembly][chrom]))

        distribution = np.random.choice(range(len(sizes)), size=nr_rand_pos, p=np.array(sizes) / np.sum(sizes))
        vals, counts = np.unique(distribution, return_counts=True)

        # then distribute inside a chromosome
        counts = [0]
        startpos = [np.array([])]
        for i, (assembly, chrom) in enumerate(combis[1:]):
            where = np.where(vals == i)
            count = where[0][0] if len(where[0]) > 0 else 0
            counts.append(count)
            startpos.append(np.random.randint(0, len(self.database.fastas[assembly][chrom]) - seq_length, size=count))

        cumsum = np.cumsum(counts)

        return combis, cumsum, startpos

    def get_onehot_sequence(self, assembly: str, chrom: str, chromstart: int, chromend: int,
                            process=None):
        """
        Get the one-hot encoded sequence based on the assembly, chromosome, chromstart and chromend.
        """
        if process:
            seq = self.fastas[process][assembly][chrom][chromstart:chromend]
        else:
            seq = self.database.fastas[assembly][chrom][chromstart:chromend]

        seq = util.sequence_to_onehot(seq)

        return seq

    @abstractmethod
    def get_label(self):
        pass

    @abstractmethod
    def get_condition(self):
        pass


class BedDataSet(DataSet):
    SELECT = " SELECT Assembly, Chromosome, ChromStart, ChromEnd "
    FROM = " FROM Bed " \
           " INNER JOIN Chromosome Chr ON Bed.ChromosomeId = Chr.ChromosomeId " \
           " INNER JOIN Condition Con  ON Bed.ConditionId  = Con.ConditionId " \
           " INNER JOIN Assembly Ass   ON Chr.AssemblyId   = Ass.AssemblyId "

    def __init__(self, database: str, where: str = "", seq_length: int = 200,
                 frac_region: float = 0.5, union_conditions: bool = True, **kwargs: int):
        DataSet.__init__(self, database, where, seq_length, **kwargs)
        self.frac_region = frac_region
        self.peaks = self.get_peak_locations(union_conditions)

    def get_condition(self):
        return None

    def get_label(self):
        pass

    def get_peak_locations(self, union) -> dict:
        """
        Map of all assemblies and genomic coordinates and whether or not it falls within a peak
        region.

        :return
        {
        assembly1:
            {chrom1: np.array([False, False, True, ..., True], dtype=bool),
             chrom2: np.array([False, True, False, ..., True], dtype=bool),
             ...
        }
        """
        peaks = dict()
        if union:
            for assembly, chrom, chromstart, chromend in self.fetchall:
                if assembly not in peaks:
                    peaks[assembly] = dict()
                if chrom not in peaks[assembly]:
                    peaks[assembly][chrom] = np.zeros(
                        len(self.database.fastas[assembly][chrom]),
                        dtype=bool
                    )
                peaks[assembly][chrom][chromstart:chromend] = True

        return peaks

    def get_label_occurence(self) -> (int, int):
        """
        Get the number of occurrences of each label.

        First returns the number of negatives (non-peak), and then the number of positives (peak).
        """
        peaks = 0
        for i, (assembly, chrom) in enumerate(self.chromosomes[1:]):
            for chromstart in self.positions[i + 1]:
                peaks += at_least(self.peaks[assembly][chrom][chromstart:chromstart+self.seq_length]
                                  , self.frac_region)

        return len(self) - peaks, peaks


# class NarrowPeakDataSet(BedDataSet, DataSet):
#     SELECT = "SELECT Assembly, Chromosome, ChromStart, Peak FROM Bed " \
#              "INNER JOIN Chromosome Chr ON Bed.ChromosomeId = Chr.ChromosomeId " \
#              "INNER JOIN Condition Con  ON Bed.ConditionId  = Con.ConditionId " \
#              "INNER JOIN Assembly Ass   ON Chr.AssemblyId   = Ass.AssemblyId "
#
#     def __init__(self, database: str, query: str = "", seq_length: int = 200, **kwargs: int):
#         raise NotImplementedError
#         DataSet.__init__(self, database, query, seq_length, **kwargs)
#         self.summits = self.get_summit_locations()
#
#     def get_summit_locations(self) -> dict:
#         """
#         Map of all assemblies and genomic coordinates and whether or not it contains a summit.
#
#         :return
#         {
#         assembly1:
#             {chrom1: np.array([False, False, True, ..., True], dtype=bool),
#              chrom2: np.array([False, True, False, ..., True], dtype=bool),
#              ...
#         }
#         """
#         summits = dict()
#         for assembly, chrom, chromstart, peak in self.fetchall:
#             if assembly not in summits:
#                 summits[assembly] = dict()
#             if chrom not in summits[assembly]:
#                 summits[assembly][chrom] = np.zeros(
#                     len(self.database.fastas[assembly][chrom]),
#                     dtype=bool
#                 )
#             summits[assembly][chrom][chromstart + peak] = True
#
#         return summits


# class BamDataSet(DataSet):
#     def __init__(self):
#         raise NotImplementedError
#
#
# class WigDataSet(DataSet):
#     def __init__(self):
#         raise NotImplementedError


import time
now = time.time()
dataset = BedDataSet("/vol/PeakSQL/peakSQL.sqlite", where="where species='danRer11'", stride=2000)


print(dataset.get_label_occurence())
print(time.time() - now)
