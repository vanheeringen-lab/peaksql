import numpy as np
from typing import List, Tuple

from .base import _DataSet


class BedDataSet(_DataSet):
    """
    The BedRegion ...
    """

    SELECT_LABEL = (
        " Bed.ConditionId, BedVirtual_{assembly}.ChromStart, "
        "BedVirtual_{assembly}.ChromEnd"
    )

    def array_from_query(
        self, query: List[Tuple[int, int, int]], chromstart: int, chromend: int,
    ) -> np.ndarray:
        positions = np.zeros((len(self.all_conditions), self.inner_range), dtype=bool)

        for condition_id, start, end in query:
            min_idx = int(start - chromstart)
            if min_idx < 0:
                min_idx = 0

            max_idx = int(end - chromstart)
            if max_idx > positions.shape[1]:
                max_idx = positions.shape[1]

            positions[condition_id, min_idx:max_idx] = True

        return positions
