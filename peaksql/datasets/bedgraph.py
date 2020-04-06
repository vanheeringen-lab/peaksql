import numpy as np
from typing import List, Tuple

from .base import _DataSet


class BedGraphDataSet(_DataSet):
    """
    The BedGraphDataSet ...
    """

    SELECT_LABEL = (
        " Bed.ConditionId, BedVirtual_{assembly}.ChromStart, "
        "BedVirtual_{assembly}.ChromEnd, Bed.DataValue "
    )

    def __init__(self, *args, **kwargs):
        kwargs.update({"label_func": "none"})
        _DataSet.__init__(self, *args, **kwargs)

    def array_from_query(
        self, query: List[Tuple[int, int, int]], chromstart: int, chromend: int,
    ) -> np.ndarray:
        positions = np.zeros(
            (len(self.all_conditions), chromend - chromstart), dtype=float
        )

        if len(query):
            condition, min_idx, max_idx, value = np.split(np.array(query), 4, axis=1)
            min_idx -= chromstart
            max_idx -= chromstart
            min_idx[min_idx < 0] = 0
            max_idx[max_idx > positions.shape[1]] = positions.shape[1]
            condition = condition.astype(int).flatten()
            min_idx = min_idx.astype(int).flatten()
            max_idx = max_idx.astype(int).flatten()
            value = value.flatten()

            positions[condition, min_idx:max_idx] = value

        return positions
