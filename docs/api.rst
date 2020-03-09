API
===============

This is a collection of public functions and classes relevant for users of PeakSQL.

DataBase
---------------

.. autoclass:: peaksql.database.DataBase
   :members: add_assembly, add_data, assemblies

DataSet loaders
---------------

peaksql.datasets.base
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: peaksql.datasets.base._DataSet
   :members: __getitem__, get_onehot_sequence, get_label

peaksql.datasets.bedregion
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: peaksql.datasets.bedregion.BedRegionDataSet
   :show-inheritance:

peaksql.datasets.narrowpeak module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: peaksql.datasets.narrowpeak.NarrowPeakDataSet
   :show-inheritance:

Util
-------------------

.. automodule:: peaksql.util
   :members: sequence_to_onehot
   :undoc-members:
   :show-inheritance:
