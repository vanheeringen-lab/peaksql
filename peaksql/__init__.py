# make sure that relevant stuff is importable
from . import database, util
from .database import DataBase
from .datasets.narrowpeak import NarrowPeakDataSet
from .datasets.bed import BedDataSet

__all__ = ["database", "util", "DataBase", "NarrowPeakDataSet", "BedDataSet"]
