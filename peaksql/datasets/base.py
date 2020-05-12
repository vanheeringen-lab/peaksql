from abc import ABC, abstractmethod
from typing import Any, Dict, List, Tuple
import math
import threading

import numpy as np
from pathos.pools import ProcessPool
import pathos

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

    # @profile
    def __init__(self, database: str, where: str = "", seq_length: int = 200, **kwargs):
        # check for valid input
        if ("stride" in kwargs) + ("nr_rand_pos" in kwargs) + (
            "nr_balanced_rand_pos" in kwargs
        ) != 1:  # xor
            raise ValueError("choose a stride OR a number of random positions")

        # store general stuff
        self.database_path = database
        self.databases: Dict[str, DataBase] = dict()
        self.seq_length = seq_length
        self.in_memory = kwargs.get("in_memory", False)
        self.iter_index = 0

        # initialize the label function
        _Labeler.__init__(self, label_func=kwargs.get("label_func", "any"))

        # sql(ite) lookup
        self.WHERE = where
        query = (
            self.SELECT_CHROM_ASS + self.FROM + self.WHERE + " ORDER BY ChromosomeId"
        )
        self._database.cursor.execute(query)
        self.fetchall = self._database.cursor.fetchall()

        # get all the conditions and their id in the database
        self.all_conditions = {
            k: v
            for k, v in self._database.cursor.execute(
                "SELECT DISTINCT Condition, ConditionId FROM Condition"
            ).fetchall()
        }

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
        if "nr_balanced_rand_pos" in kwargs:
            self.nr_balanced_rand_pos = kwargs["nr_balanced_rand_pos"]
            (
                self.chromosomes,
                self.cumsum,
                self.positions,
            ) = self.get_balanced_random_positions(
                self.seq_length, self.nr_balanced_rand_pos
            )

        # mark fetchall for garbage collection (large and we don't need it anymore)
        del self.fetchall

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
        process = pathos.helpers.mp.current_process().name + str(
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
        chromstart) triplet. The positions are sampled across the query with an even
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

    # @profile
    def get_random_positions(
        self, seq_length: int, nr_rand_pos: int
    ) -> Tuple[list, np.ndarray, np.ndarray]:
        """
        Calculate a map that connects __getitem__ indices to (assembly, chrom,
        chromstart) triplet. The positions are sampled across the query randomly, but
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

    def get_balanced_random_positions(
        self, seq_length: int, nr_rand_pos: int, cores: int = None
    ) -> Tuple[list, np.ndarray, np.ndarray]:
        """
        Calculate a map that connects __getitem__ indices to (assembly, chrom,
        chromstart) triplet. The positions are sampled across the query randomly, but
        proportional to the size of each chromosome.

        The first return value is a list of (assembly, chrom) pairs, the second list
        consists of the cumulative sum of the number of indices that belong to this
        assembly, chrom pair, and the third list contains all chromstarts of the
        sequences. This allows for a decently fast and memory-efficient lookup of
        genomic positions corresponding to an index.
        """
        get_label = getattr(self, "get_label", None)
        if not callable(get_label):
            raise NotImplementedError("Can not label ")

        if cores is None:
            cores = pathos.helpers.mp.cpu_count()

        # setup a processing pool
        pool = ProcessPool(nodes=cores)

        true_goal = false_goal = nr_rand_pos // 2 + (nr_rand_pos % 2 > 0)
        d_non_empty_combis: Dict[Tuple[Any, ...], List] = {(None, None): []}
        trues = falses = total = 1

        while trues + falses < nr_rand_pos:
            # we overestimate a bit for faster convergence
            nr_pos = int((nr_rand_pos * (total / trues) - total) * 1.05)
            if nr_pos > 1_000_000:
                nr_pos = 1_000_000

            s_non_empty_combis, _, s_startpos = self.get_random_positions(
                seq_length, nr_pos
            )
            s_endpos = [start + seq_length for start in s_startpos]

            args = []
            for (assembly, chrom), chromstarts, chromends in zip(
                s_non_empty_combis, s_startpos, s_endpos
            ):
                for chromstart, chromend in zip(chromstarts, chromends):
                    args.append((assembly, chrom, chromstart, chromend))

            # turns args into a list of chunks to spread across cores
            chunksize = math.ceil(nr_pos / cores)
            chunk_args = [
                args[i : i + chunksize] for i in range(0, len(args), chunksize)
            ]
            with util.MinimizedDataset(self):
                labels = list(pool.map(self.__labels_from_array, chunk_args))

            for i, arg in enumerate(args):
                label = labels[i // chunksize][i % chunksize]
                if label[0] and trues <= true_goal:
                    d_non_empty_combis.setdefault(tuple(arg[:2]), []).append(arg[2])
                    trues += 1
                elif label[0] == False and falses <= false_goal:  # noqa: E712
                    d_non_empty_combis.setdefault(tuple(arg[:2]), []).append(arg[2])
                    falses += 1
                total += 1

        startpos = list(d_non_empty_combis.values())
        cumsum = np.cumsum([len(positions) for positions in startpos])
        non_empty_combis = list(d_non_empty_combis.keys())
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
                  ({chromend} >= BedVirtual_{assembly}.ChromStart)
        """.format(
            assembly=assembly
        )
        query_result = self._database.cursor.execute(query).fetchall()
        positions = self.array_from_query(query_result, chromstart, chromend)
        labels = self.label_from_array(positions)

        return labels

    def downsample(self, ncpus: int = None):
        """

        """
        raise NotImplementedError
        if ncpus is None:
            ncpus = 1

        trues = []
        falses = []
        for i in range(len(self)):
            site = self._index_to_site(i)
            label = self.get_label(*site)
            if label[0]:
                trues.append(i)
            else:
                falses.append(i)

        lowest = min([len(trues), len(falses)])
        trues = np.random.choice(trues, lowest, replace=False)
        falses = np.random.choice(falses, lowest, replace=False)

        indices = np.sort(np.concatenate((trues, falses)))
        new_chromosomes = [(None, None)]
        new_positions = [[]]
        ass_chrom_idx = 1
        for index in indices:
            while index >= self.cumsum[ass_chrom_idx]:
                ass_chrom_idx += 1

            if self.chromosomes[ass_chrom_idx] not in new_chromosomes:
                new_chromosomes.append(self.chromosomes[ass_chrom_idx])
                new_positions.append([])

            new_positions[-1].append(
                self.positions[ass_chrom_idx][index - self.cumsum[ass_chrom_idx - 1]]
            )

        self.chromosomes = new_chromosomes
        self.positions = new_positions
        self.cumsum = np.cumsum([len(pos) for pos in new_positions])

    def upsample(self):
        """

        """
        raise NotImplementedError

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

    def __labels_from_array(self, positions: np.ndarray) -> np.ndarray:
        """

        """
        return [self.get_label(*pos) for pos in positions]
