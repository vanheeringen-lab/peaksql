import numpy as np
from typing import List, Tuple

from .base import _DataSet


class BedGraphDataSet(_DataSet):
    """
    The BedGraphDataSet ...
    """

    SELECT_LABEL = (
        " Bed.ConditionId, BedVirtual.ChromStart, Bed.ChromEnd, Bed.DataValue "
    )

    def __init__(self, **kwargs):
        _DataSet.__init__(self, **kwargs)
        self.label_from_array = self.none

    def array_from_query(
        self, query: List[Tuple[int, int, int]], chromstart: int, chromend: int,
    ) -> np.ndarray:
        positions = np.zeros(
            (len(self.all_conditions), chromend - chromstart), dtype=float
        )

        for condition_id, start, end, value in query:
            min_idx = int(start - chromstart)
            if min_idx < 0:
                min_idx = 0

            max_idx = int(end - chromstart)
            if max_idx > positions.shape[1]:
                max_idx = positions.shape[1]

            positions[condition_id, min_idx:max_idx] = value

        return positions
