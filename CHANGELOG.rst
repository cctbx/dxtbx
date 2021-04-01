DIALS 3.4.1 (2021-04-01)
========================

Bugfixes
--------

- Nexus: Diamond Light Source beamlines are now properly identified (`#339 <https://github.com/cctbx/dxtbx/issues/339>`_)


DIALS 3.4.0 (2021-03-15)
========================

Features
--------

- ``FormatHDF5SaclaMPCCD`` is now a "Lazy load" format (`#227 <https://github.com/cctbx/dxtbx/issues/227>`_)
- Show image counts when displaying ``Scan`` objects (e.g. ``dials.show``) (`#271 <https://github.com/cctbx/dxtbx/issues/271>`_)
- The ``Scan.append`` default tolerance is increased to 3% of the image width, to
  accommodate electron diffraction datasets with poor rotation stages. (`#277 <https://github.com/cctbx/dxtbx/issues/277>`_)
- Preliminary support for images derived from Timepix 2M detector in NeXus / NXmx format (`#298 <https://github.com/cctbx/dxtbx/issues/298>`_)
- Add function ``dxtbx.util.get_url_scheme``, to identify URL-style image paths in a cross-platform way (`#301 <https://github.com/cctbx/dxtbx/issues/301>`_)
- Add support for raw data from the SwissFEL Jungfrau 16M detector, including some estimates of pixel errors (`#303 <https://github.com/cctbx/dxtbx/issues/303>`_)
- CBF decompression: Validate expected image size, and the ``cbf_decompress``
  function now accepts the output array size, and returns the number of
  items read. (`#313 <https://github.com/cctbx/dxtbx/issues/313>`_)
- Include test for equality of ``PxMmStrategy`` in ``Panel`` equality operator. (`#319 <https://github.com/cctbx/dxtbx/issues/319>`_)
- Format suport for Eiger 16M XE at Diamond - recognise legacy and updated beamline names. (`#323 <https://github.com/cctbx/dxtbx/issues/323>`_)
- The function ``ExperimentList.from_templates`` has been added for construction convenience (`#333 <https://github.com/cctbx/dxtbx/issues/333>`_)


Bugfixes
--------

- Fix Gatan DM4 format reader. (`#297 <https://github.com/cctbx/dxtbx/issues/297>`_)
- Fix ``dxtbx.`` commands crashing on Windows when unicode output is directed to a file (`#306 <https://github.com/cctbx/dxtbx/issues/306>`_)
- ``dxtbx.dlsnxs2cbf``: Properly display help message when passed ``-h`` (`#309 <https://github.com/cctbx/dxtbx/issues/309>`_)
- Check for existence of certain numpy types before using them. (`#318 <https://github.com/cctbx/dxtbx/issues/318>`_)
- Correctly link to HDF5 shared libraries on Windows (`#329 <https://github.com/cctbx/dxtbx/issues/329>`_)


Deprecations and Removals
-------------------------

- The main development branch of dxtbx was renamed from 'master' to 'main'. (`#281 <https://github.com/cctbx/dxtbx/issues/281>`_)
- ``DataBlock`` is now deprecated. Please use ``ExperimentList`` instead. (`#288 <https://github.com/cctbx/dxtbx/issues/288>`_)
- Remove obsolete format ``FormatNexusExternalDataFile`` (`#328 <https://github.com/cctbx/dxtbx/issues/328>`_)
- The previously deprecated ``ScanFactory.single`` has been removed. Use ``ScanFactory.single_file`` instead. (`#332 <https://github.com/cctbx/dxtbx/issues/332>`_)
- ``ExperimentListTemplateImporter`` is now deprecated. Please use ``ExperimentList.from_templates``. (`#333 <https://github.com/cctbx/dxtbx/issues/333>`_)


Misc
----

- `#272 <https://github.com/cctbx/dxtbx/issues/272>`_, `#275 <https://github.com/cctbx/dxtbx/issues/275>`_, `#279 <https://github.com/cctbx/dxtbx/issues/279>`_, `#282 <https://github.com/cctbx/dxtbx/issues/282>`_, `#287 <https://github.com/cctbx/dxtbx/issues/287>`_, `#288 <https://github.com/cctbx/dxtbx/issues/288>`_, `#291 <https://github.com/cctbx/dxtbx/issues/291>`_, `#293 <https://github.com/cctbx/dxtbx/issues/293>`_, `#302 <https://github.com/cctbx/dxtbx/issues/302>`_, `#308 <https://github.com/cctbx/dxtbx/issues/308>`_, `#316 <https://github.com/cctbx/dxtbx/issues/316>`_, `#320 <https://github.com/cctbx/dxtbx/issues/320>`_, `#322 <https://github.com/cctbx/dxtbx/issues/322>`_, `#324 <https://github.com/cctbx/dxtbx/issues/324>`_, `#326 <https://github.com/cctbx/dxtbx/issues/326>`_, `#327 <https://github.com/cctbx/dxtbx/issues/327>`_, `#331 <https://github.com/cctbx/dxtbx/issues/331>`_


DIALS 3.3.4 (2021-03-05)
========================

Bugfixes
--------

- Fix error corrupting data when writing CBF files with large pixel values.
  This affected ``dxtbx.dlsnxs2cbf`` and ``dials.merge_cbf`` (`#314 <https://github.com/cctbx/dxtbx/issues/314>`_)


DIALS 3.3.3 (2021-02-15)
========================

Bugfixes
--------

- Fix for missing ``SENSOR_THICKNESS=`` in XDS.INP generated for EIGER datasets introduced in 3.3.1 (`#296 <https://github.com/cctbx/dxtbx/issues/296>`_)


DIALS 3.3.2 (2021-02-01)
========================

Bugfixes
--------

- Don't interpret windows paths as URIs, causing failure to import images (`#284 <https://github.com/cctbx/dxtbx/issues/284>`_)
- Fix bug in ``nexus.DataFactory`` that allowed access to twice as many
  images as available on disk for VDS nexus files. (`#285 <https://github.com/cctbx/dxtbx/issues/285>`_)
- Bug fix for live per-image analysis of HDF5/SWMR files, ensuring that
  a process can see data for images written after a process first sees
  a given data file. (`#289 <https://github.com/cctbx/dxtbx/issues/289>`_)
- Bug fix for generating XDS.INP for eiger datasets - ensure that
  ``DETECTOR=EIGER (not PILATUS)`` (`#292 <https://github.com/cctbx/dxtbx/issues/292>`_)


DIALS 3.3.1 (2021-01-18)
========================

Features
--------

- NeXus files are now opened in SWMR mode. (`#270 <https://github.com/cctbx/dxtbx/issues/270>`_)


DIALS 3.3.0 (2021-01-04)
========================

Features
--------

- ``FormatMultiImage``: When constructing an imageset with the indices of some
  (not all) single images in the container, we skip reading models for the
  images that were not requested. In some cases this speeds up imageset
  construction by 8x. (`#210 <https://github.com/cctbx/dxtbx/issues/210>`_)
- Read detector distance from the XTC streams for LCLS Jungfrau data (`#246 <https://github.com/cctbx/dxtbx/issues/246>`_)
- Set the per-shot gain for the ePix and Jungfrau detectors at LCLS. (`#250 <https://github.com/cctbx/dxtbx/issues/250>`_)
- Allow format classes to be marked as ``@abstract``. This means that they will
  be considered and returned by the Registry search if they are the best match,
  but are intended to represent an incomplete "category" of format class that
  other classes build on, so cannot be instantiated. (`#255 <https://github.com/cctbx/dxtbx/issues/255>`_)


Bugfixes
--------

- When creating "Lazy" ImageSets the static mask from the image file was not being properly applied (`#227 <https://github.com/cctbx/dxtbx/issues/227>`_)
- Be more robust when handling nexus scan axes (`#252 <https://github.com/cctbx/dxtbx/issues/252>`_)
- Improve error message when attempting to import data-only h5 files (`#261 <https://github.com/cctbx/dxtbx/issues/261>`_)
- Fix finding HDF5 plugins when using dials-installer (`#265 <https://github.com/cctbx/dxtbx/issues/265>`_)
- Prevent errors reading eiger data, if ``h5py`` is imported before dxtbx (`#266 <https://github.com/cctbx/dxtbx/issues/266>`_)
- Fix errors introduced by moving to ``h5py`` 3.1+ (`#267 <https://github.com/cctbx/dxtbx/issues/267>`_)
- Improve error message when attempting to import unsupported files (`#1220 <https://github.com/cctbx/dxtbx/issues/1220>`_)


Deprecations and Removals
-------------------------

- Deprecate ``ScanFactory.single``. Please use ``ScanFactory.single_file``
  without the `format=` argument, which has been removed. `ScanFactory.single`
  will be removed in a future version. (`#233 <https://github.com/cctbx/dxtbx/issues/233>`_)
- Remove deprecated ``dxtbx.serialize.dump.experiment_list``, ``dxtbx.serialize.filename.load_path``,
  and ``as_str`` argument to ``dxtbx.serialize.xds.to_xds().XDS_INP()`` (`#248 <https://github.com/cctbx/dxtbx/issues/248>`_)
- The ``ignore()`` functionality on Format classes has been removed. Such
  classes should be marked as ``@abstract`` instead. (`#255 <https://github.com/cctbx/dxtbx/issues/255>`_)
- Deprecate the HDF5 plugin discovery patch that is applied when dxtbx is
  imported before h5py. Please update your HDF5 plugins package. (`#258 <https://github.com/cctbx/dxtbx/issues/258>`_)
- Remove ``FormatHDF5RawData`` format class. This was only ever used
  experimentally, and caused confusion when incorrectly importing nexus
  side files. (`#261 <https://github.com/cctbx/dxtbx/issues/261>`_)
- The deprecated ``dxtbx.datablock.DataBlockDumper`` and ``serialize.dump``
  have been removed. (`#269 <https://github.com/cctbx/dxtbx/issues/269>`_)


Misc
----

- `#238 <https://github.com/cctbx/dxtbx/issues/238>`_, `#257 <https://github.com/cctbx/dxtbx/issues/257>`_, `#260 <https://github.com/cctbx/dxtbx/issues/260>`_, `#262 <https://github.com/cctbx/dxtbx/issues/262>`_, `#267 <https://github.com/cctbx/dxtbx/issues/267>`_


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
