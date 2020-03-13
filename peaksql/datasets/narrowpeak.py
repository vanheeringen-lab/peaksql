import numpy as np
from typing import List, Tuple

from .labeler import _Labeler
from .base import _DataSet


class NarrowPeakDataSet(_DataSet, _Labeler):
    """
    The NarrowPeakDataSet ...
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
        _DataSet.__init__(self, database, where, seq_length, **kwargs)
        _Labeler.__init__(self, label_func=label_func)

    def array_from_query(
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
