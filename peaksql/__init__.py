__version__ = '0.0.0'

# make sure that relevant stuff is importable
from . import database, util
from .database import DataBase
from .datasets.narrowpeak import NarrowPeakDataSet
from .datasets.bedregion import BedRegionDataSet
