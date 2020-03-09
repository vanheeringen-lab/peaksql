import numpy as np
from abc import ABC, abstractmethod
from typing import List, Tuple

from .base import _DataSet


class _BedDataSet(_DataSet, ABC):
    """
    The BedDataSet...
    """

    FROM = (
        " FROM Chromosome Chr "
        " INNER JOIN Assembly Ass  ON Chr.AssemblyId   = Ass.AssemblyId "
    )
    SELECT_LABEL: str

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

        setattr(self, "_label_from_array", eval("self.label_" + kwargs["label_func"]))

    def label_any(self, positions: np.ndarray) -> np.ndarray:
        return np.any(positions, axis=1)

    def label_inner_any(self, positions: np.ndarray) -> np.ndarray:
        mid = positions.shape[0] // 2
        return self.label_any(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    def label_all(self, positions: np.ndarray) -> np.ndarray:
        return np.all(positions, axis=1)

    def label_inner_all(self, positions: np.ndarray) -> np.ndarray:
        mid = positions.shape[0] // 2
        return self.label_all(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    def label_fraction(self, positions: np.ndarray) -> np.ndarray:
        return np.sum(positions, axis=1) / positions.shape[1] >= self.fraction

    def label_inner_fraction(self, positions: np.ndarray) -> np.ndarray:
        mid = positions.shape[0] // 2
        return self.label_fraction(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    @abstractmethod
    def _array_from_query(
        self,
        query: List[Tuple[int, int, int, int]],
        cur_chrom_id: int,
        chromstart: int,
        chromend: int,
    ) -> np.ndarray:
        pass

    def _label_from_array(self, positions: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def get_label(
        self, assembly: str, chrom: str, chromstart: int, chromend: int
    ) -> np.ndarray:
        """
        Get the label that corresponds to chromstart:chromend.
        """
        assemblyid = self._database.get_assembly_id(assembly)
        chromosomeid = self._database.get_chrom_id(assemblyid, chrom)

        offset = self._database.cursor.execute(
            f"""
            SELECT Offset FROM Chromosome WHERE ChromosomeId = {chromosomeid}
            """
        ).fetchone()[0]
        chromstart += offset
        chromend += offset

        query = f"""
            SELECT {self.SELECT_LABEL}
            FROM BedVirtual
            INNER JOIN Bed on BedVirtual.BedId = Bed.BedId
            WHERE ({chromstart} <= BedVirtual.ChromEnd) AND
                  ({chromend} >= BedVirtual.ChromStart)
        """
        query_result = self._database.cursor.execute(query).fetchall()
        positions = self._array_from_query(
            query_result, chromosomeid, chromstart, chromend
        )
        labels = self._label_from_array(positions)

        return labels
