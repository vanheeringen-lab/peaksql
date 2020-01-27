# Should database be able to generate these? db.get_NarrowPeakDataSet(query)
# answer: probably yes, TODO

# TODO: make work for different conditions
import numpy as np
import multiprocessing

import pyfaidx

from .util import sequence_to_onehot, _binary_search, jit_any


class NarrowPeakDataSet:
    """
    Dataset ...
    """
    SELECT = "SELECT Assembly, Chromosome, ChromStart, ChromStop FROM Bed " \
             "INNER JOIN Chromosome Chr ON Bed.ChromosomeId = Chr.ChromosomeId " \
             "INNER JOIN Condition Con  ON Bed.ConditionId  = Con.ConditionId " \
             "INNER JOIN Assembly Ass   ON Chr.AssemblyId   = Ass.AssemblyId "

    def __init__(self, database, query, seq_length, stride):
        """

        """
        # store ...
        self.seq_length = seq_length
        self.stride = stride
        self.database = database
        self.fastas = {}

        query = self.SELECT + query
        database.cursor.execute(query)
        self.fetchall = database.cursor.fetchall()

        # TODO proper test
        assert len(self.fetchall[0]) == 4, "TODO"

        self.peaks = self.get_peak_locations()
        self.order, self.lens = self.get_size_map()

    def __len__(self) -> int:
        """

        """
        return self.lens[-1]

    def __getitem__(self, index: int):
        """
        TODO: type hint, docstring
        """
        process = self._get_process()
        assembly, chrom, chromstart, chromstop = self._index_to_site(index)

        # get the sequence
        seq = self.fastas[process][assembly][chrom][chromstart:chromstop]
        seq = sequence_to_onehot(seq)

        # get the label
        label = np.any(self.peaks[assembly][chrom][chromstart:chromstop])

        return seq, label

    def _index_to_site(self, index: int) -> (str, str, int, int):
        """
        Convert the index of self.__getitem__ to a tuple of (assembly, chrom, chromstart, chromstop)

        Uses binary search for fast retrieval.
        """
        index_bs = _binary_search(index, self.lens)
        pos = (index - self.lens[index_bs - 1]) * self.stride
        assembly, chrom = self.order[index_bs]

        return assembly, chrom, pos, pos + self.seq_length

    def _get_process(self):
        """
        PyFaidx is not multiprocessing safe when reading from fasta index. However if we start a
        Fasta class for each process, the file pointers don't get mangled, and we can use
        multiprocessing. This allows us to "stream" our data in e.g. a Pytorch dataloader.
        """
        process = multiprocessing.current_process().name
        if process not in self.fastas:
            self.database.cursor.execute("SELECT Assembly, AbsPath FROM Assembly")
            self.fastas[process] = {
                assembly: pyfaidx.Fasta(abspath)
                for assembly, abspath in self.database.cursor.fetchall()
            }
        return process

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

    def get_size_map(self) -> (list, list):
        """
        Build a map that connects __getitem__ indices to (assembly, chrom) pairs.

        The first return value is a list of (assembly, chrom) pairs, and the second list consists of
        the cumulative sum of the number of indices that belong to this pair. The second list can
        then be used to quickly find a genomic position that belongs to an index.
        """
        order = [(None, None)]
        lens = [0]

        for assembly in self.peaks:
            for chrom in self.peaks[assembly]:
                lens.append((len(self.peaks[assembly][chrom]) - self.seq_length) // self.stride)
                order.append((assembly, chrom))
        lens = np.cumsum(lens)
        return order, lens

    def get_label_occurence(self) -> (int, int):
        """
        Get the number of occurrences of each label.

        First returns the number of
        """
        peaks = 0
        for i in range(len(self)):
            assembly, chrom, chromstart, chromstop = self._index_to_site(i)
            peaks += jit_any(self.peaks[assembly][chrom][chromstart:chromstop])

        return len(self) - peaks, peaks
