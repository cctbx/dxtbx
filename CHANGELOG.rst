DIALS 2.0 (2019-10-23)
======================

Features
--------

- Change dxtbx format registry to using entry points

  dxtbx now discovers format classes during configuration time instead of
  at runtime. Format classes can either be added into the dxtbx/format
  directory as usual, registered by other python packages using the
  'dxtbx.format' entry point, or installed by the user via the
  'dxtbx.install_format' command.

  To register format classes stored in ~/.dxtbx you need to run
  'dxtbx.install_format -u' whenever you add or remove format classes.

  Changes for library users:
  * A number of registry lookup methods were deprecated or removed.
  * Exceptions from format .understand() methods are no longer discarded.
    Similarly, when no matching format was found the datablock find_format()
    methods now return 'None' and no longer raise exceptions.
    In both cases the caller will need to deal with the situation appropriately.
  * Format classes must be named 'Format*', and must inherit either from
    other format classes or from the top-level format class, 'Format'.
    Base classes must be given as their original name and must therefore not
    contain '.'s. (#34)
- Reading compressed FullCBF files - .gz or .bz2 - is now supported (#72)
- Add an optional Format.get_static_mask() method

  This allows format classes to define a static mask to be used across all images
  in an imageset. (#73)
- Add new command dxtbx.dlsnxs2cbf which converts Nexus files created at
  Diamond Light Source to .cbf files. (#81)
- Added ``ExperimentList.from_file`` for easily loading data. This means
  that experiment lists and reflection tables can now load the same way. (#100)


Bugfixes
--------

- Replace h5py `visititems` with `local_visit` implementation to work around using soft links in Eiger / hdf5 files. (#75)
- Fix FormatNexusEigerDLS16M.understand() for 2019/run4 datasets (#85)
- Reduce number of redundant file operations in dxtbx

  This includes a change in the DataBlock() construction semantics: sequences from
  identical detectors are merged into a single DataBlock() object regardless of
  their position in the call order. Since DataBlock() is deprecated and any
  reliance on order would have to be handled explicitly downstream anyway this
  should not have any impact on users or developers. (#89)
- Fix setting a per-panel pedestal

  Per-panel pedestals are now respected when the corrected data is used. (#108)


Misc
----

- #76, #90


Dxtbx 2.0 (2019-10-23)
======================

Features
--------

- Change dxtbx format registry to using entry points

  dxtbx now discovers format classes during configuration time instead of
  at runtime. Format classes can either be added into the dxtbx/format
  directory as usual, registered by other python packages using the
  'dxtbx.format' entry point, or installed by the user via the
  'dxtbx.install_format' command.

  To register format classes stored in ~/.dxtbx you need to run
  'dxtbx.install_format -u' whenever you add or remove format classes.

  Changes for library users:
  * A number of registry lookup methods were deprecated or removed.
  * Exceptions from format .understand() methods are no longer discarded.
    Similarly, when no matching format was found the datablock find_format()
    methods now return 'None' and no longer raise exceptions.
    In both cases the caller will need to deal with the situation appropriately.
  * Format classes must be named 'Format*', and must inherit either from
    other format classes or from the top-level format class, 'Format'.
    Base classes must be given as their original name and must therefore not
    contain '.'s. (#34)
- Reading compressed FullCBF files - .gz or .bz2 - is now supported (#72)
- Add an optional Format.get_static_mask() method

  This allows format classes to define a static mask to be used across all images
  in an imageset. (#73)
- Add new command dxtbx.dlsnxs2cbf which converts Nexus files created at
  Diamond Light Source to .cbf files. (#81)
- Added ``ExperimentList.from_file`` for easily loading data. This means
  that experiment lists and reflection tables can now load the same way. (#100)


Bugfixes
--------

- Replace h5py `visititems` with `local_visit` implementation to work around using soft links in Eiger / hdf5 files. (#75)
- Fix FormatNexusEigerDLS16M.understand() for 2019/run4 datasets (#85)
- Reduce number of redundant file operations in dxtbx

  This includes a change in the DataBlock() construction semantics: sequences from
  identical detectors are merged into a single DataBlock() object regardless of
  their position in the call order. Since DataBlock() is deprecated and any
  reliance on order would have to be handled explicitly downstream anyway this
  should not have any impact on users or developers. (#89)
- Fix setting a per-panel pedestal

  Per-panel pedestals are now respected when the corrected data is used. (#108)


Misc
----

- #76, #90
