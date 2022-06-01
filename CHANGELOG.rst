DIALS 3.8.5 (2022-06-01)
========================

Features
--------

- When installed as a libtbx module, dxtbx will not install python packages into ``conda_base/``. (`#511 <https://github.com/cctbx/dxtbx/issues/511>`_)


Bugfixes
--------

- ``FormatCBFFullPilatus``: Handle detector information better in cases of multiple or missing panels. (`#508 <https://github.com/cctbx/dxtbx/issues/508>`_)
- Remove check for beam/normalization orthogonality in ``Beam.rotate_around_origin``. This could stop processing of older, incorrectly configured data. (`#510 <https://github.com/cctbx/dxtbx/issues/510>`_)
- Correct a unicode error reading Bruker ``.sfrm`` files. With thanks to `Dennis Brookner <https://github.com/dennisbrookner>`_ for this change. (`#518 <https://github.com/cctbx/dxtbx/issues/518>`_)


Misc
----

- `#513 <https://github.com/cctbx/dxtbx/issues/513>`_, `#520 <https://github.com/cctbx/dxtbx/issues/520>`_


DIALS 3.8.4 (2022-04-01)
========================

Bugfixes
--------

- ``FormatNXmxI19_2``:  Allow data from beamline I19-2 at Diamond Light Source to be processed with optional masking of the beamline's standard diamond anvil pressure cell with a 76° aperture. (`#481 <https://github.com/cctbx/dxtbx/issues/481>`_)


DIALS 3.8.3 (2022-02-22)
========================

Bugfixes
--------

- FormatNXmx: Open nexus files in SWMR mode. (`#478 <https://github.com/cctbx/dxtbx/issues/478>`_)


DIALS 3.8.2 (2022-02-07)
========================

Bugfixes
--------

- ``dxtbx.dlsnxs2cbf``: Provide more general support for correctly formatted NXmx-flavoured NeXus data.  Previously, only a very limited subset of experiment geometries and data formats were supported. (`#453 <https://github.com/cctbx/dxtbx/issues/453>`_)
- More robustly handle different ways of recording single-value NXmx detector metadata. (`#460 <https://github.com/cctbx/dxtbx/issues/460>`_)
- Fix ``dxtbx.plot_detector_models`` running on newer matplotlib versions. (`#475 <https://github.com/cctbx/dxtbx/issues/475>`_)


DIALS 3.8.1 (2022-01-25)
========================

Features
--------

- Updated bad pixel mask for DLS I23 PILATUS 12M for 2022 run 1 (`#469 <https://github.com/cctbx/dxtbx/issues/469>`_)


dxtbx 3.8.0 (2022-01-11)
========================

Features
--------

- ``dxtbx`` can be optionally used without ``cbflib_adaptbx``. (`#368 <https://github.com/cctbx/dxtbx/issues/368>`_)
- Experimental support for building dxtbx with CMake. (`#449 <https://github.com/cctbx/dxtbx/issues/449>`_)
- Track dxtbx version explicitly, with bump2version. (`#458 <https://github.com/cctbx/dxtbx/issues/458>`_)


Bugfixes
--------

- Fix an arithmetic mistake in ``dxtbx.model.Goniometer.rotate_around_origin``, which was mangling the addition of a new rotation to the goniostat rotation operator :math:`\mathbf{R}`. (`#451 <https://github.com/cctbx/dxtbx/issues/451>`_)
- Correct pedestal handling for simulated images from ``simtbx`` (`#456 <https://github.com/cctbx/dxtbx/issues/456>`_)
- Ensure ``FormatTIFF`` only understands images with the expected basic TIFF header. (`#457 <https://github.com/cctbx/dxtbx/issues/457>`_)
- Get CI builds working again by restricting setuptools<60. (`#459 <https://github.com/cctbx/dxtbx/issues/459>`_)


Improved Documentation
----------------------

- Update the documentation of the in-house convention for representing the goniostat rotation operator :math:`\mathbf{R}`, to match `the conventions page <https://dials.github.io/documentation/conventions.html#the-dxtbx-goniometer-model>`_ of the online DIALS documentation. (`#450 <https://github.com/cctbx/dxtbx/issues/450>`_)


Deprecations and Removals
-------------------------

- Remove ``ImageToEwaldSphere``, which was used in a now-removed utility. (`#446 <https://github.com/cctbx/dxtbx/issues/446>`_)
- The deprecated function ``dxtbx.model.detector_helpers.project_2d`` has been removed. The deprecation warning on usage of `DataBlock` has been made more visible. (`#448 <https://github.com/cctbx/dxtbx/issues/448>`_)


Misc
----

- `#366 <https://github.com/cctbx/dxtbx/issues/366>`_


DIALS 3.7.0 (2021-11-01)
========================

Features
--------

- New function ``Crystal.clone()``, to get a new Crystal object of the same type. (`#420 <https://github.com/cctbx/dxtbx/issues/420>`_)
- New ``fast_slow_beam_centre=`` parameter for detector models allows setting the beam centre using fast, slow [panel] value ordering. (`#421 <https://github.com/cctbx/dxtbx/issues/421>`_)
- Added ``dlstbx.nexus.nxmx`` module that provides a high-level read-only interface to HDF5 files adhering to the NeXus/NXmx standard, and support for Diamond Light Source's I19-2 EIGER detector. (`#423 <https://github.com/cctbx/dxtbx/issues/423>`_)
- Allow importing experiment lists from single-file templates. (`#425 <https://github.com/cctbx/dxtbx/issues/425>`_)
- Support NeXus data from the Tristan event-mode detector on beamline I19 at Diamond Light Source. (`#428 <https://github.com/cctbx/dxtbx/issues/428>`_)


Bugfixes
--------

- Fix installation using Python 3.7 on Windows. (`#441 <https://github.com/cctbx/dxtbx/issues/441>`_)
- Better support for detector SMV ADSC SN442. (`#445 <https://github.com/cctbx/dxtbx/issues/445>`_)


Deprecations and Removals
-------------------------

- The function ``dxtbx.model.detector_helpers.project_2d`` has been renamed ``get_detector_projection_2d_axes``. Usage of the function ``project_2d`` is deprecated and will be removed after DIALS 3.7. (`#422 <https://github.com/cctbx/dxtbx/issues/422>`_)
- Drop support for Python 3.6. (`#424 <https://github.com/cctbx/dxtbx/issues/424>`_)


Misc
----

- `#394 <https://github.com/cctbx/dxtbx/issues/394>`_, `#422 <https://github.com/cctbx/dxtbx/issues/422>`_, `#430 <https://github.com/cctbx/dxtbx/issues/430>`_, `#431 <https://github.com/cctbx/dxtbx/issues/431>`_, `#432 <https://github.com/cctbx/dxtbx/issues/432>`_, `#435 <https://github.com/cctbx/dxtbx/issues/435>`_, `#436 <https://github.com/cctbx/dxtbx/issues/436>`_


DIALS 3.6.2 (2021-09-21)
========================

Bugfixes
--------

- Fix broken ``dxtbx.install_format`` command. (`#434 <https://github.com/cctbx/dxtbx/issues/434>`_)


DIALS 3.6.0 (2021-08-16)
========================

Features
--------

- Add **experimental** ``dxtbx.flumpy.to_numpy``, ``.from_numpy``, ``.vec_from_numpy`` and
  ``.mat3_from_numpy`` for zero-copy conversions between numpy and `scitbx.array_family.flex``
  arrays. There is also a lower-level class ``Scuffer`` that allows exposing of flex arrays via
  generic python buffer interfaces for e.g. Cython interoperability. (`#377 <https://github.com/cctbx/dxtbx/issues/377>`_)
- ``ExperimentListFactory.from_filenames(...)``, ``Format.get_imageset(...)``, and
  ``ImageSetFactory.new(...)`` now accept objects implementing the Python file system path protocol
  (PEP-519). (`#386 <https://github.com/cctbx/dxtbx/issues/386>`_)


Bugfixes
--------

- Fix support of older FormatSMVADSCSN442 images (`#369 <https://github.com/cctbx/dxtbx/issues/369>`_)
- More detailed error messages are now printed after internal ``H5Dread`` calls fail (`#374 <https://github.com/cctbx/dxtbx/issues/374>`_)
- Fix error reading BioMAX data with H5py 3.3 (`#389 <https://github.com/cctbx/dxtbx/issues/389>`_)
- Fix potential problem where mask geometry was unfixable (`#411 <https://github.com/cctbx/dxtbx/issues/411>`_)
- Handle installing dxtbx as a "real" package when the ``conda_base/`` is read-only (`#413 <https://github.com/cctbx/dxtbx/issues/413>`_)
- Check for empty beams in XTC streams (`#419 <https://github.com/cctbx/dxtbx/issues/419>`_)


Deprecations and Removals
-------------------------

- The previously deprecated ``ExperimentListTemplateImporter`` has been removed. Please use
  ``ExperimentList.from_templates`` instead. (`#333 <https://github.com/cctbx/dxtbx/issues/333>`_)


Misc
----

- Move dxtbx to ``src/`` layout, and install as a package (`#382 <https://github.com/cctbx/dxtbx/pull/382>`_)
- `#311 <https://github.com/cctbx/dxtbx/issues/311>`_, `#373 <https://github.com/cctbx/dxtbx/issues/373>`_, `#375 <https://github.com/cctbx/dxtbx/issues/375>`_, `#380 <https://github.com/cctbx/dxtbx/issues/380>`_, `#381 <https://github.com/cctbx/dxtbx/issues/381>`_, `#384 <https://github.com/cctbx/dxtbx/issues/384>`_, `#386 <https://github.com/cctbx/dxtbx/issues/386>`_, `#388 <https://github.com/cctbx/dxtbx/issues/388>`_, `#390 <https://github.com/cctbx/dxtbx/issues/390>`_, `#391 <https://github.com/cctbx/dxtbx/issues/391>`_, `#396 <https://github.com/cctbx/dxtbx/issues/396>`_, `#400 <https://github.com/cctbx/dxtbx/issues/400>`_, `#401 <https://github.com/cctbx/dxtbx/issues/401>`_, `#402 <https://github.com/cctbx/dxtbx/issues/402>`_, `#403 <https://github.com/cctbx/dxtbx/issues/403>`_, `#404 <https://github.com/cctbx/dxtbx/issues/404>`_


DIALS 3.5.4 (2021-07-27)
========================

Bugfixes
--------

- Allow reading of new SACLA hdf5 data (`#408 <https://github.com/cctbx/dxtbx/issues/408>`_)


DIALS 3.5.2 (2021-06-28)
========================

Bugfixes
--------

- End the I03 "bad mask" duration, since it is now masked at the file level. (`#385 <https://github.com/cctbx/dxtbx/issues/385>`_)
- ``dxtbx.dlsnxs2cbf``: Handle missing chi/phi axis entries. (`#387 <https://github.com/cctbx/dxtbx/issues/387>`_)


DIALS 3.5.1 (2021-06-14)
========================

Bugfixes
--------

- Extend duration of bad module mask for Diamond I03 EIGER 2XE 16M detector indefinitely. This will be updated in a future release. (`#370 <https://github.com/cctbx/dxtbx/issues/370>`_)
- Handle scan data which wraps through 0° instead of >=360° (`#379 <https://github.com/cctbx/dxtbx/issues/379>`_)


DIALS 3.5.0 (2021-05-27)
========================

Features
--------

- Add ``FormatMRC.py`` for electron diffraction images and image stacks recorded on Thermo Fisher microscopes (`#335 <https://github.com/cctbx/dxtbx/issues/335>`_)
- Improved support for Gatan DM4 format images and stacks (`#338 <https://github.com/cctbx/dxtbx/issues/338>`_)
- Improved support for TIA (Emispec) .ser files (`#345 <https://github.com/cctbx/dxtbx/issues/345>`_)
- Improved support for ``.emi`` sidecar files in ``FormatSER`` (`#354 <https://github.com/cctbx/dxtbx/issues/354>`_)
- Add support for Python 3.9. (`#365 <https://github.com/cctbx/dxtbx/issues/365>`_)


Bugfixes
--------

- Bug fixes for extended header reading in ``FormatMRC.py`` (`#343 <https://github.com/cctbx/dxtbx/issues/343>`_)
- ``dxtbx.dlsnxs2cbf``: Fixed on Windows using ``hdf5plugin`` (`#344 <https://github.com/cctbx/dxtbx/issues/344>`_)
- Mask temporarily bad modules on the Diamond I03 EIGER 2XE 16M detector (`#348 <https://github.com/cctbx/dxtbx/issues/348>`_)
- Fix rare error during CBF compression (`#352 <https://github.com/cctbx/dxtbx/issues/352>`_)
- Extend duration of bad module mask for Diamond I03 EIGER 2XE 16M detector (`#355 <https://github.com/cctbx/dxtbx/issues/355>`_)


Deprecations and Removals
-------------------------

- Remove legacy HDF5 plugin handling. Please update your conda environment if you still have issues. (`#340 <https://github.com/cctbx/dxtbx/issues/340>`_)
- Remove classes and functions deprecated in the previous release: ``dxtbx.datablock.*Diff``, ``dxtbx.model.experiment_list.SequenceDiff``, ``dxtbx.serialize.load.imageset_from_string``. (`#347 <https://github.com/cctbx/dxtbx/issues/347>`_)
- Removed unused support for reading experiments from pickle files (`#361 <https://github.com/cctbx/dxtbx/issues/361>`_)
- Remove the ability to save experiments in pickle format (`#363 <https://github.com/cctbx/dxtbx/issues/363>`_)


Misc
----

- `#334 <https://github.com/cctbx/dxtbx/issues/334>`_, `#337 <https://github.com/cctbx/dxtbx/issues/337>`_, `#342 <https://github.com/cctbx/dxtbx/issues/342>`_, `#346 <https://github.com/cctbx/dxtbx/issues/346>`_, `#350 <https://github.com/cctbx/dxtbx/issues/350>`_, `#351 <https://github.com/cctbx/dxtbx/issues/351>`_, `#353 <https://github.com/cctbx/dxtbx/issues/353>`_, `#357 <https://github.com/cctbx/dxtbx/issues/357>`_, `#360 <https://github.com/cctbx/dxtbx/issues/360>`_, `#364 <https://github.com/cctbx/dxtbx/issues/364>`_


DIALS 3.4.1 (2021-03-31)
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
- Format support for Eiger 16M XE at Diamond - recognise legacy and updated beamline names. (`#323 <https://github.com/cctbx/dxtbx/issues/323>`_)
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
