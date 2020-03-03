import numpy as np
import multiprocessing
from abc import ABC, abstractmethod
from typing import Dict, Tuple

from ..database import DataBase
import peaksql.util as util


class _DataSet(ABC):
    """
    DataSet baseclass.
    """

    SELECT_CHROM_ASS = "SELECT Assembly, Chromosome "
    FROM: str

    def __init__(
        self, database: str, where: str = "", seq_length: int = 200, **kwargs: int
    ):
        # check for valid input
        if ("stride" in kwargs) == ("nr_rand_pos" in kwargs):  # xor
            raise ValueError("choose a stride OR a number of random positions")

        # store general stuff
        self.database_path = database
        self.databases: Dict[str, DataBase] = dict()
        self.seq_length = seq_length
        self.in_memory = kwargs.get("in_memory", False)

        # sql(ite) lookup
        self.WHERE = where
        query = (
            self.SELECT_CHROM_ASS + self.FROM + self.WHERE + " ORDER BY ChromosomeId"
        )
        self.database.cursor.execute(query)
        self.fetchall = self.database.cursor.fetchall()

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
            for k, v in self.database.cursor.execute(
                "SELECT DISTINCT Condition, ConditionId FROM Condition"
            ).fetchall()
        }

        # mark fetchall for garbage collection (large and we don't need it anymore)
        del self.fetchall

    def __len__(self) -> int:
        """
        Return the number of indices this dataset contains.
        """
        return self.cumsum[-1]

    def __getitem__(self, index: int) -> Tuple[np.array, int]:
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

    def _get_process(self):
        """
        PyFaidx is not multiprocessing safe when reading from fasta index or with
        sql(ite) queries. However if we start a new DataBase class for each process, we
        automatically start new sql(ite) connections and pyfaidx instances. This allows
        us to "stream" our data parallel in e.g. a Pytorch dataloader.
        """
        process = multiprocessing.current_process().name
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

    def get_strided_positions(self, seq_length: int, stride: int):
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
                0, len(self.database.fastas[assembly][chrom]) - seq_length + 1, stride
            )
            if len(positions):
                non_empty_combis.append((assembly, chrom))
                startpos.append(positions)
                counts.append(len(startpos[-1]))

        cumsum = np.cumsum(counts)

        return non_empty_combis, cumsum, startpos

    def get_random_positions(self, seq_length, nr_rand_pos):
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
        # FIXME: this implementation is very buggy, let's just raise error for now
        combis = list({(assembly, chrom) for assembly, chrom, *_ in self.fetchall})
        combis = [(None, None)] + [
            (assembly, chrom)
            for assembly, chrom, *_ in self.fetchall
            if len(self.database.fastas[assembly][chrom]) > seq_length
        ]

        # distribute the positions over the chromosomes
        sizes = []
        for assembly, chrom in combis[1:]:
            sizes.append(len(self.database.fastas[assembly][chrom]))

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
                        len(self.database.fastas[assembly][chrom]) - seq_length,
                        size=counts[where[0][0]],
                    )
                )
                non_empty_combis.append((assembly, chrom))

        cumsum = np.cumsum(total_counts)

        return non_empty_combis, cumsum, startpos

    def get_onehot_sequence(
        self, assembly: str, chrom: str, chromstart: int, chromend: int
    ):
        """
        Get the one-hot encoded sequence based on the assembly, chromosome, chromstart
        and chromend.
        """
        process = self._get_process()

        seq = self.databases[process].fastas[assembly][chrom][chromstart:chromend]
        seq = util.sequence_to_onehot(seq)

        return seq

    @abstractmethod
    def get_label(self, *args):
        pass

    @property
    def database(self):
        process = self._get_process()
        return self.databases[process]


class _BedDataSet(_DataSet, ABC):
    """
    The BedDataSet...
    """

    FROM = (
        " FROM Chromosome Chr "
        " INNER JOIN Assembly Ass  ON Chr.AssemblyId   = Ass.AssemblyId "
    )

    def __init__(self, database: str, where: str = "", seq_length: int = 200, **kwargs):
        _DataSet.__init__(self, database, where, seq_length, **kwargs)

        assert "label_func" in kwargs and kwargs["label_func"] in [
            "any",
            "inner_any",
            "all",
            "inner_all",
            "fraction",
            "inner_fraction",
        ]

        if "inner" in kwargs["label_func"]:
            if "inner_range" not in kwargs:
                raise ValueError(
                    f"You specified an 'inner' function {kwargs['label_func']} but did "
                    f"not specify an inner_range."
                )
            self.inner_range = kwargs["inner_range"]

        if "fraction" in kwargs["label_func"]:
            if "fraction" not in kwargs:
                raise ValueError(
                    f"You specified a 'fraction' of the sequence to be within a ragion,"
                    f" but you did not not specify the fraction."
                )
            self.fraction = kwargs["fraction"]

        setattr(self, "label_from_array", eval("self.label_" + kwargs["label_func"]))

    def label_any(self, positions):
        return np.any(positions, axis=1)

    def label_inner_any(self, positions):
        mid = positions.shape[0] // 2
        return self.label_any(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    def label_all(self, positions):
        return np.all(positions, axis=1)

    def label_inner_all(self, positions):
        mid = positions.shape[0] // 2
        return self.label_all(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    def label_fraction(self, positions):
        return np.sum(positions, axis=1) / positions.shape[1] >= self.fraction

    def label_inner_fraction(self, positions):
        mid = positions.shape[0] // 2
        return self.label_fraction(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    @abstractmethod
    def array_from_query(self):
        pass

    def get_label(self, assembly, chrom, chromstart, chromend):
        """
        Get the label that corresponds to chromstart:chromend.
        """
        assemblyid = self.database.get_assembly_id(assembly)
        chromosomeid = self.database.get_chrom_id(assemblyid, chrom)

        bed_virtual = f"BedVirtual_{assemblyid}"
        query = f"""
            SELECT {self.SELECT_LABEL} FROM {bed_virtual}
            INNER JOIN Bed on {bed_virtual}.BedId = Bed.BedId
            WHERE ({chromstart} <= {bed_virtual}.ChromEnd) AND
                  ({chromend} >= {bed_virtual}.ChromStart)
        """
        query_result = self.database.cursor.execute(query).fetchall()

        positions = self.array_from_query(
            query_result, chromosomeid, chromstart, chromend
        )
        labels = self.label_from_array(positions)

        return labels
