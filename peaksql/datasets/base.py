import numpy as np
import multiprocessing
import threading
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

from ..database import DataBase
from .labeler import _Labeler
import peaksql.util as util


class _DataSet(ABC, _Labeler):
    """
    DataSet baseclass.
    """

    SELECT_CHROM_ASS = "SELECT Assembly, Chromosome "
    FROM = (
        " FROM Chromosome Chr "
        " INNER JOIN Assembly Ass  ON Chr.AssemblyId   = Ass.AssemblyId "
    )
    SELECT_LABEL: str

    def __init__(self, database: str, where: str = "", seq_length: int = 200, **kwargs):
        # check for valid input
        if ("stride" in kwargs) == ("nr_rand_pos" in kwargs):  # xor
            raise ValueError("choose a stride OR a number of random positions")

        # store general stuff
        self.database_path = database
        self.databases: Dict[str, DataBase] = dict()
        self.seq_length = seq_length
        self.in_memory = kwargs.get("in_memory", False)
        self.iter_index = 0

        # sql(ite) lookup
        self.WHERE = where
        query = (
            self.SELECT_CHROM_ASS + self.FROM + self.WHERE + " ORDER BY ChromosomeId"
        )
        self._database.cursor.execute(query)
        self.fetchall = self._database.cursor.fetchall()

        # get the genomic positions of our indices
        if "stride" in kwargs:
            self.stride = kwargs["stride"]
            self.chromosomes, self.cumsum, self.positions = self.get_strided_positions(
                self.seq_length, self.stride
            )
        if "nr_rand_pos" in kwargs:
            self.nr_rand_pos = kwargs["nr_rand_pos"]
            self.chromosomes, self.cumsum, self.positions = self.get_random_positions(
                self.seq_length, self.nr_rand_pos
            )

        # get all the conditions and their id in the database
        self.all_conditions = {
            k: v
            for k, v in self._database.cursor.execute(
                "SELECT DISTINCT Condition, ConditionId FROM Condition"
            ).fetchall()
        }

        # mark fetchall for garbage collection (large and we don't need it anymore)
        del self.fetchall

        _Labeler.__init__(self, label_func=kwargs.get("label_func", "any"))

    def __len__(self) -> int:
        """
        Return the number of indices this dataset contains.
        """
        return self.cumsum[-1]

    def __getitem__(self, index: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return the sequence in one-hot encoding and the label of the corresponding
        index.
        """
        if index >= len(self):
            raise StopIteration

        assembly, chrom, chromstart, chromend = self._index_to_site(index)

        # get the sequence, label and condition
        seq = self.get_onehot_sequence(assembly, chrom, chromstart, chromend)
        label = self.get_label(assembly, chrom, chromstart, chromend)

        return seq, label

    def _get_process(self) -> str:
        """
        PyFaidx is not multiprocessing safe when reading from fasta index or with
        sql(ite) queries. However if we start a new DataBase class for each process, we
        automatically start new sql(ite) connections and pyfaidx instances. This allows
        us to "stream" our data parallel in e.g. a Pytorch dataloader.
        """
        process = multiprocessing.current_process().name + str(
            threading.current_thread().ident
        )
        if process not in self.databases:
            self.databases[process] = DataBase(
                self.database_path, in_memory=self.in_memory
            )
        return process

    def _index_to_site(self, index: int) -> Tuple[str, str, int, int]:
        """
        Convert the index of self.__getitem__ to a tuple of (assembly, chrom,
        chromstart, chromend)

        Uses binary search for fast retrieval.
        """
        index_bs = util.binary_search(index, self.cumsum)

        assembly, chrom = self.chromosomes[index_bs]
        chromstart = self.positions[index_bs][index - self.cumsum[index_bs - 1]]
        chromend = chromstart + self.seq_length

        return assembly, chrom, chromstart, chromend

    def get_strided_positions(
        self, seq_length: int, stride: int
    ) -> Tuple[list, np.ndarray, np.ndarray]:
        """
        Calculate a map that connects __getitem__ indices to (assembly, chrom,
        chromstart) triplet. The positions are sampled accross the query with an even
        stride.

        The first return value is a list of (assembly, chrom) pairs, the second list
        consists of the cumulative sum of the number of indices that belong to this
        assembly, chrom pair, and the third list contains all chromstarts of the
        sequences. This allows for a decently fast and memory-efficient lookup of
        genomic positions corresponding to an index.
        """
        combis = [(None, None)] + list(
            {(assembly, chrom) for assembly, chrom, *_ in self.fetchall}
        )

        counts = [0]
        startpos = [np.array([])]
        non_empty_combis = [(None, None)]
        for assembly, chrom in combis[1:]:
            positions = np.arange(
                0, len(self._database.fastas[assembly][chrom]) - seq_length + 1, stride
            )
            if len(positions):
                non_empty_combis.append((assembly, chrom))
                startpos.append(positions)
                counts.append(len(startpos[-1]))

        cumsum = np.cumsum(counts)

        return non_empty_combis, cumsum, startpos

    def get_random_positions(
        self, seq_length: int, nr_rand_pos: int
    ) -> Tuple[list, np.ndarray, np.ndarray]:
        """
        Calculate a map that connects __getitem__ indices to (assembly, chrom,
        chromstart) triplet. The positions are sampled accross the query randomly, but
        proportional to the size of each chromosome.

        The first return value is a list of (assembly, chrom) pairs, the second list
        consists of the cumulative sum of the number of indices that belong to this
        assembly, chrom pair, and the third list contains all chromstarts of the
        sequences. This allows for a decently fast and memory-efficient lookup of
        genomic positions corresponding to an index.
        """
        combis = [(None, None)] + [
            (assembly, chrom)
            for assembly, chrom, *_ in self.fetchall
            if len(self._database.fastas[assembly][chrom]) > seq_length
        ]

        # distribute the positions over the chromosomes
        sizes = []
        for assembly, chrom in combis[1:]:
            sizes.append(len(self._database.fastas[assembly][chrom]))

        distribution = np.random.choice(
            range(len(sizes)), size=nr_rand_pos, p=np.array(sizes) / np.sum(sizes)
        )
        vals, counts = np.unique(distribution, return_counts=True)

        # then distribute inside a chromosome
        total_counts = [0]
        startpos = [np.array([])]
        non_empty_combis = [(None, None)]
        for i, (assembly, chrom) in enumerate(combis[1:]):
            where = np.where(vals == i)
            if len(where[0]) > 0:
                total_counts.append(counts[where[0][0]])
                startpos.append(
                    np.random.randint(
                        0,
                        len(self._database.fastas[assembly][chrom]) - seq_length,
                        size=counts[where[0][0]],
                    )
                )
                non_empty_combis.append((assembly, chrom))

        cumsum = np.cumsum(total_counts)

        return non_empty_combis, cumsum, startpos

    def get_onehot_sequence(
        self, assembly: str, chrom: str, chromstart: int, chromend: int
    ) -> np.ndarray:
        """
        Get the one-hot encoded sequence based on the assembly, chromosome, chromstart
        and chromend.
        """
        seq = self._database.fastas[assembly][chrom][chromstart:chromend]
        seq = util.sequence_to_onehot(seq)

        return seq

    def get_label(
        self, assembly: str, chrom: str, chromstart: int, chromend: int
    ) -> np.ndarray:
        """
        Get the label that corresponds to chromstart:chromend.
        """
        offset, chromosomeid = self._database.get_offset_chromosomeid(assembly, chrom)
        chromstart += offset
        chromend += offset

        # if we only want the label of an inner part
        if hasattr(self, "inner_range"):
            midpoint = chromstart + self.seq_length // 2
            chromstart = midpoint - self.inner_range
            chromend = midpoint + self.inner_range

        query = f"""
            SELECT {self.SELECT_LABEL}
            FROM BedVirtual_{assembly}
            INNER JOIN Bed on BedVirtual_{assembly}.BedId = Bed.BedId
            WHERE ({chromstart} < BedVirtual_{assembly}.ChromEnd) AND
                  ({chromend} >= BedVirtual_{assembly}.ChromStart) AND
                  ChromosomeId = {chromosomeid}
        """.format(
            assembly=assembly
        )
        query_result = self._database.cursor.execute(query).fetchall()
        positions = self.array_from_query(query_result, chromstart, chromend)
        labels = self.label_from_array(positions)

        return labels

    @property
    def _database(self):
        process = self._get_process()
        return self.databases[process]

    @abstractmethod
    def array_from_query(
        self, query: List[Tuple[int, int, int]], chromstart: int, chromend: int,
    ) -> np.ndarray:
        pass

    def label_from_array(self, positions: np.ndarray) -> np.ndarray:
        raise NotImplementedError
