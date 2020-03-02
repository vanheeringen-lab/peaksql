import numpy as np

from .basedataset import _BedDataSet


class NarrowPeakDataSet(_BedDataSet):
    """
    The NarrowPeakDataSet ...
    """

    SELECT_LABEL = " Bed.ChromosomeId, Bed.ConditionId, Bed.ChromStart, Bed.Peak "

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

        for chromosome_id, condition_id, start, peak in query:
            if cur_chrom_id == chromosome_id:
                peak_idx = start - chromstart + peak
                if 0 <= peak_idx < positions.shape[1]:
                    positions[condition_id, peak_idx] = True

        return positions
