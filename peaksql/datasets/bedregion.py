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
                if condition_id is None:
                    condition_idx = 0
                else:
                    condition_idx = self.all_conditions.index(condition_id)
                min_idx = start - chromstart if start - chromstart >= 0 else 0
                max_idx = (
                    end - chromstart
                    if end - chromstart <= positions.shape[1]
                    else positions.shape[1]
                )
                positions[condition_idx, min_idx:max_idx] = True

        return positions
