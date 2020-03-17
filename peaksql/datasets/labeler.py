import numpy as np


class _Labeler:
    """
    Labeler...
    """

    def __init__(self, **kwargs):

        print(kwargs["label_func"])
        assert "label_func" in kwargs and kwargs["label_func"] in [
            "any",
            "all",
            "fraction",
        ]

        self.inner_range = kwargs.get("inner_range", self.seq_length)
        self.ratio = kwargs.get("ratio", 1.0)

        setattr(self, "label_from_array", eval("self." + kwargs["label_func"]))

    def any(self, positions: np.ndarray) -> np.ndarray:
        return np.any(positions, axis=1)

    def all(self, positions: np.ndarray) -> np.ndarray:
        return np.all(positions, axis=1)

    def fraction(self, positions: np.ndarray) -> np.ndarray:
        return np.sum(positions, axis=1) / positions.shape[1] >= self.ratio
