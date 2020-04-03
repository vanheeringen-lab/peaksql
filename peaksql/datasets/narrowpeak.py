import numpy as np
from typing import List, Tuple

from .base import _DataSet


class NarrowPeakDataSet(_DataSet):
    """
    The NarrowPeakDataSet expects that narrowPeak files have been added to the DataBase.
    """

    SELECT_LABEL = (
        " Bed.ChromosomeId, Bed.ConditionId, BedVirtual_{assembly}.ChromStart, Bed.Peak"
    )

    def array_from_query(
        self, query: List[Tuple[int, int, int]], chromstart: int, chromend: int,
    ) -> np.ndarray:
        positions = np.zeros((len(self.all_conditions), self.inner_range), dtype=bool)

        for condition_id, start, peak in query:
            peak_idx = int(start - chromstart + peak)
            if 0 <= peak_idx < positions.shape[1]:
                positions[condition_id, peak_idx] = True

        return positions
