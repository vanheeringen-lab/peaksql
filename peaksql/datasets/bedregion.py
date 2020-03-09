import numpy as np
from abc import ABC
from typing import List, Tuple

from .bed import _BedDataSet


class BedRegionDataSet(_BedDataSet, ABC):
    """
    The BedRegion DataSet
    """

    SELECT_LABEL = (
        " Bed.ChromosomeId, Bed.ConditionId, BedVirtual.ChromStart, BedVirtual.ChromEnd"
    )

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

        for chromosome_id, condition_id, start, end in query:
            if cur_chrom_id == chromosome_id:
                min_idx = max(0, int(start - chromstart))
                max_idx = min(positions.shape[1], int(end - chromstart))

                positions[condition_id, min_idx:max_idx] = True

        return positions
