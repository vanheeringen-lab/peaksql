import numpy as np

from .basedataset import _BedDataSet


class BedRegionDataSet(_BedDataSet):
    """
    The BedRegion ...
    """

    SELECT_LABEL = " Bed.ChromosomeId, Bed.ConditionId, Bed.ChromStart, Bed.ChromEnd "

    def __init__(
        self,
        database: str,
        where: str = "",
        seq_length: int = 200,
        label_func: str = "any",
        **kwargs: int
    ):
        _BedDataSet.__init__(
            self, database, where, seq_length, label_func=label_func, **kwargs
        )

    def array_from_query(self, query, cur_chrom_id, chromstart, chromend):
        positions = np.zeros(
            (len(self.all_conditions), chromend - chromstart), dtype=bool
        )

        for chromosome_id, condition_id, start, end in query:
            if cur_chrom_id == chromosome_id:
                min_idx = max(0, start - chromstart)
                max_idx = min(positions.shape[1], end - chromstart)

                positions[condition_id, min_idx:max_idx] = True

        return positions
