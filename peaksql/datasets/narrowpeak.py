import numpy as np
from abc import ABC
from typing import List, Tuple

from .bed import _BedDataSet


class NarrowPeakDataSet(_BedDataSet, ABC):
    """
    The NarrowPeakDataSet expects that narrowPeak files have been added to the DataBase.
    """

    SELECT_LABEL = " Bed.ChromosomeId, Bed.ConditionId, BedVirtual.ChromStart, Bed.Peak"

    def __init__(
        self,
        database: str,
        where: str = "",
        seq_length: int = 200,
        label_func: str = "any",
        **kwargs
    ):
        _BedDataSet.__init__(
            self, database, where, seq_length, label_func=label_func, **kwargs
        )

    def _array_from_query(
        self,
        query: List[Tuple[int, int, int, int]],
        cur_chrom_id: int,
        chromstart: int,
        chromend: int,
    ) -> np.ndarray:
        positions = np.zeros(
            (len(self.all_conditions), chromend - chromstart), dtype=bool
        )

        for chromosome_id, condition_id, start, peak in query:
            if cur_chrom_id == chromosome_id:
                peak_idx = start - chromstart + peak
                if 0 <= peak_idx < positions.shape[1]:
                    positions[condition_id, peak_idx] = True

        return positions
