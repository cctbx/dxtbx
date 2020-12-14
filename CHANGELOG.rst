DIALS 3.2.0 (2020-10-27)
========================

Features
--------

- Add ``clear_cache()`` method to clear internal imageset cache  (`#218 <https://github.com/cctbx/dxtbx/issues/218>`_)
- Add ``dxtbx.model.detector_helpers.project_2d`` function, which calculates
  a 2D projection of the detector panels into a frame aligned to the
  image. This is intended for use in display tasks for mostly co-planar
  detectors.  (`#224 <https://github.com/cctbx/dxtbx/issues/224>`_)
- image template: add support for ``nameNNNN`` e.g. ``image1234`` as a valid name  (`#234 <https://github.com/cctbx/dxtbx/issues/234>`_)
- ``BeamFactory.simple`` will now return an unpolarised beam for >~247 KeV beams
  (e.g. Electron diffraction)  (`#243 <https://github.com/cctbx/dxtbx/issues/243>`_)


Bugfixes
--------

- Fix reading of legacy pickle-image files created from Python 3  (`#205 <https://github.com/cctbx/dxtbx/issues/205>`_)
- Allow importing filenames with special format characters like ``%``  (`#214 <https://github.com/cctbx/dxtbx/issues/214>`_)
- ``dxtbx.dlsnxs2cbf``: strip timezone when making CBF file timestamps  (`#235 <https://github.com/cctbx/dxtbx/issues/235>`_)
- Fix error reading nexus files when using hardlinks to detector models  (`#240 <https://github.com/cctbx/dxtbx/issues/240>`_)
- SMV Formats: Use header gain values if present, rather than guessing  (`#242 <https://github.com/cctbx/dxtbx/issues/242>`_)


Misc
----
- We have moved the pytest launchers from cctbx_project to dxtbx. If you run
  into ``libtbx.configure`` errors make sure both repositories are up to date  (`#231 <https://github.com/cctbx/dxtbx/issues/231>`_)
- `#209 <https://github.com/cctbx/dxtbx/issues/209>`_, `#211 <https://github.com/cctbx/dxtbx/issues/211>`_,
  `#212 <https://github.com/cctbx/dxtbx/issues/212>`_, `#217 <https://github.com/cctbx/dxtbx/issues/217>`_,
  `#225 <https://github.com/cctbx/dxtbx/issues/225>`_, `#226 <https://github.com/cctbx/dxtbx/issues/226>`_,
  `#230 <https://github.com/cctbx/dxtbx/issues/230>`_


DIALS 3.1.4 (2020-10-12)
========================

Bugfixes
--------

- Handle more errors using Eiger-Nexus files


DIALS 3.1.3 (2020-09-28)
========================

Bugfixes
--------

- ``dxtbx.image_average``: Better use of MPI to avoid errors and increase
  performance  (`#207 <https://github.com/cctbx/dxtbx/issues/207>`_)
- Update DLS I23 bad pixel mask after detector has been cleaned, fixing
  previously bad modules.  (`#220 <https://github.com/cctbx/dxtbx/issues/220>`_)
- Change default bit depth for DLS eigers where header information is missing


DIALS 3.1.1 (2020-09-01)
========================

Bugfixes
--------

- Don't crash handling FormatSMVADSC images with floating-point pedestal values  (`#216 <https://github.com/cctbx/dxtbx/issues/216>`_)
- Allow importing filenames with special format characters like %  (`#214 <https://github.com/cctbx/dxtbx/issues/214>`_)


DIALS 3.1 (2020-08-17)
======================

Features
--------

- Add generic multi-panel support for FormatCBFMiniPilatus and subclasses. Data
  matching format classes inheriting from FormatCBFMiniPilatus can now be
  imported with the option multi_panel=True to treat the detector as multiple
  panels, instead of a single panel comprising the whole detector.  (`#177 <https://github.com/cctbx/dxtbx/issues/177>`_)
- New tool ``dxtbx.show_mask_info`` to show the number of masked pixels for each module  (`#198 <https://github.com/cctbx/dxtbx/issues/198>`_)
- **Experimental - Alpha API**: Add Spectrum as a read-only class obtainable from
  an imageset, and implement reading spectra from NeXus files.  (`#201 <https://github.com/cctbx/dxtbx/issues/201>`_)


Bugfixes
--------

- Better handle string conversion when NeXus files  (`#190 <https://github.com/cctbx/dxtbx/issues/190>`_)
- HDF5 / NeXus: Correctly use the mask if available.  (`#198 <https://github.com/cctbx/dxtbx/issues/198>`_)


DIALS 3.0.4 (2020-07-20)
========================

- HDF5 / NeXus: Read image dimensions directly from dataset shape instead of
  reported image_size, as latter can sometimes be backwards  (`#189 <https://github.com/cctbx/dxtbx/issues/189>`_)
- Support image_range when importing images into an ImageSet so only a subset
  of the images are used
- Diamond-specific Eiger/Nexus: Fix handling of masked pixels in the image so
  that module join regions are no longer marked as overloaded (i.e. yellow) in
  the image viewer  (`#180 <https://github.com/cctbx/dxtbx/issues/180>`_)


DIALS 3.0.2 (2020-06-23)
========================

Bugfixes
--------

- Fix sensor-material handling for Jungfrau 4M and 16M detectors


DIALS 3.0.1 (2020-06-11)
========================

Bugfixes
--------

- Account for beam centre record changing with ADSC 442 move from 8.3.1 to 5.0.1  (`#171 <https://github.com/cctbx/dxtbx/issues/171>`_)
- Fix handling for hierarchical NeXus detectors  (`#175 <https://github.com/cctbx/dxtbx/issues/175>`_)
- Prevent mangling of URL-based filenames via abspath  (`#176 <https://github.com/cctbx/dxtbx/issues/176>`_)
- Fix incorrect axis detection on MAX IV Eiger and Spring8  (`#178 <https://github.com/cctbx/dxtbx/issues/178>`_)


DIALS 3.0 (2020-05-18)
======================

Features
--------

- A new recalculated unit cell attribute is added to the Crystal model, for use by post-integration cell refinement methods, such as that of dials.two_theta_refine.  (`#142 <https://github.com/cctbx/dxtbx/issues/142>`_)
- Add ExperimentList.change_basis() convenience method.  (`#166 <https://github.com/cctbx/dxtbx/issues/166>`_)
- Allow creation of Format classes that accept URLs instead of files  (`#173 <https://github.com/cctbx/dxtbx/issues/173>`_)


Bugfixes
--------

- Fix a bug whereby reading a single-image data set from an Eiger detector would lead to an error.  (`#156 <https://github.com/cctbx/dxtbx/issues/156>`_)
- Fix formatting of unit cell parameters with negligible standard uncertainties  (`#165 <https://github.com/cctbx/dxtbx/issues/165>`_)
- New Eiger FileWriter (20.1.16.56035) produces NeXus compliant files, which exposed a bug in finding axis sample depends on, now fixed.  (`#168 <https://github.com/cctbx/dxtbx/issues/168>`_)


Misc
----

- `#164 <https://github.com/cctbx/dxtbx/issues/164>`_


DIALS 2.2 (2020-03-15)
======================

Bugfixes
--------

- Fix spot-finding on images with file names ending in '0000.cbf'  (`#133 <https://github.com/cctbx/dxtbx/issues/133>`_)
- Fixed imageset slicing for image sets starting from image 0  (`#141 <https://github.com/cctbx/dxtbx/issues/141>`_)


DIALS 2.1 (2019-12-16)
======================

Features
--------

- With changes in dials.import sequences of stills are imported as individual
  experiments all dereferencing one image set - this is the change set to support
  this on load.  (`#118 <https://github.com/cctbx/dxtbx/issues/118>`_)


Bugfixes
--------

- Reinstate support for historic VMXi EIGER 1 images  (`#119 <https://github.com/cctbx/dxtbx/issues/119>`_)
- Fix crash when opening dataset containing many .h5 files  (`#126 <https://github.com/cctbx/dxtbx/issues/126>`_)


Deprecations and Removals
-------------------------

- dxtbx extensions can no longer be imported from `dxtbx`
  and must now be imported from `dxtbx.ext`  (`#29 <https://github.com/cctbx/dxtbx/issues/29>`_)


Misc
----

- `#124 <https://github.com/cctbx/dxtbx/issues/124>`_


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
    contain '.'s.  (`#34 <https://github.com/cctbx/dxtbx/issues/34>`_)
- Reading compressed FullCBF files - .gz or .bz2 - is now supported  (`#72 <https://github.com/cctbx/dxtbx/issues/72>`_)
- Add an optional Format.get_static_mask() method

  This allows format classes to define a static mask to be used across all images
  in an imageset.  (`#73 <https://github.com/cctbx/dxtbx/issues/73>`_)
- Add new command dxtbx.dlsnxs2cbf which converts Nexus files created at
  Diamond Light Source to .cbf files.  (`#81 <https://github.com/cctbx/dxtbx/issues/81>`_)
- Added ``ExperimentList.from_file`` for easily loading data. This means
  that experiment lists and reflection tables can now load the same way.  (`#100 <https://github.com/cctbx/dxtbx/issues/100>`_)


Bugfixes
--------

- Replace h5py `visititems` with `local_visit` implementation to work around using soft links in Eiger / hdf5 files.  (`#75 <https://github.com/cctbx/dxtbx/issues/75>`_)
- Fix FormatNexusEigerDLS16M.understand() for 2019/run4 datasets  (`#85 <https://github.com/cctbx/dxtbx/issues/85>`_)
- Reduce number of redundant file operations in dxtbx

  This includes a change in the DataBlock() construction semantics: sequences from
  identical detectors are merged into a single DataBlock() object regardless of
  their position in the call order. Since DataBlock() is deprecated and any
  reliance on order would have to be handled explicitly downstream anyway this
  should not have any impact on users or developers.  (`#89 <https://github.com/cctbx/dxtbx/issues/89>`_)
- Fix setting a per-panel pedestal

  Per-panel pedestals are now respected when the corrected data is used.  (`#108 <https://github.com/cctbx/dxtbx/issues/108>`_)


Misc
----

- `#76 <https://github.com/cctbx/dxtbx/issues/76>`_, `#90 <https://github.com/cctbx/dxtbx/issues/90>`_
