# make sure that relevant stuff is importable
from . import database, util
from .database import DataBase
from .datasets.narrowpeak import NarrowPeakDataSet
from .datasets.bedregion import BedRegionDataSet

__all__ = ["database", "util", "DataBase", "NarrowPeakDataSet", "BedRegionDataSet"]
