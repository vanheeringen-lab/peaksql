import numpy as np


class _Labeler:
    """
    Labeler...
    """

    FROM = (
        " FROM Chromosome Chr "
        " INNER JOIN Assembly Ass  ON Chr.AssemblyId   = Ass.AssemblyId "
    )
    SELECT_LABEL: str

    def __init__(self, **kwargs):

        assert "label_func" in kwargs and kwargs["label_func"] in [
            "any",
            "inner_any",
            "all",
            "inner_all",
            "fraction",
            "inner_fraction",
        ]

        if "inner" in kwargs["label_func"]:
            if "inner_range" not in kwargs:
                raise ValueError(
                    f"You specified an 'inner' function {kwargs['label_func']} but did "
                    f"not specify an inner_range."
                )
            self.inner_range = kwargs["inner_range"]

        if "fraction" in kwargs["label_func"]:
            if "fraction" not in kwargs:
                raise ValueError(
                    f"You specified a 'fraction' of the sequence to be within a region,"
                    f" but you did not not specify the fraction."
                )
            self.fraction = kwargs["fraction"]

        setattr(self, "label_from_array", eval("self.label_" + kwargs["label_func"]))

    def label_any(self, positions: np.ndarray) -> np.ndarray:
        return np.any(positions, axis=1)

    def label_inner_any(self, positions: np.ndarray) -> np.ndarray:
        mid = positions.shape[0] // 2
        return self.label_any(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    def label_all(self, positions: np.ndarray) -> np.ndarray:
        return np.all(positions, axis=1)

    def label_inner_all(self, positions: np.ndarray) -> np.ndarray:
        mid = positions.shape[0] // 2
        return self.label_all(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )

    def label_fraction(self, positions: np.ndarray) -> np.ndarray:
        return np.sum(positions, axis=1) / positions.shape[1] >= self.fraction

    def label_inner_fraction(self, positions: np.ndarray) -> np.ndarray:
        mid = positions.shape[0] // 2
        return self.label_fraction(
            positions[:, mid - self.inner_range : mid + self.inner_range + 1]
        )
