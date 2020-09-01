DIALS 3.1.1 (2020-09-01)
========================

Bugfixes
--------

- Don't crash handling FormatSMVADSC images with floating-point pedestal values (#216)
- Allow importing filenames with special format characters like % (#214)


DIALS 3.1 (2020-08-17)
======================

Features
--------

- Add generic multi-panel support for FormatCBFMiniPilatus and subclasses. Data
  matching format classes inheriting from FormatCBFMiniPilatus can now be
  imported with the option multi_panel=True to treat the detector as multiple
  panels, instead of a single panel comprising the whole detector. (#177)
- New tool ``dxtbx.show_mask_info`` to show the number of masked pixels for each module (#198)
- **Experimental - Alpha API**: Add Spectrum as a read-only class obtainable from
  an imageset, and implement reading spectra from NeXus files. (#201)


Bugfixes
--------

- Better handle string conversion when NeXus files (#190)
- HDF5 / NeXus: Correctly use the mask if available. (#198)


DIALS 3.0.4 (2020-07-20)
========================

- HDF5 / NeXus: Read image dimensions directly from dataset shape instead of
  reported image_size, as latter can sometimes be backwards (#189)
- Support image_range when importing images into an ImageSet so only a subset
  of the images are used
- Diamond-specific Eiger/Nexus: Fix handling of masked pixels in the image so
  that module join regions are no longer marked as overloaded (i.e. yellow) in
  the image viewer (#180)


DIALS 3.0.2 (2020-06-23)
========================

Bugfixes
--------

- Fix sensor-material handling for Jungfrau 4M and 16M detectors


DIALS 3.0.1 (2020-06-11)
========================

Bugfixes
--------

- Account for beam centre record changing with ADSC 442 move from 8.3.1 to 5.0.1 (#171)
- Fix handling for hierarchical NeXus detectors (#175)
- Prevent mangling of URL-based filenames via abspath (#176)
- Fix incorrect axis detection on MAX IV Eiger and Spring8 (#178)


DIALS 3.0 (2020-05-18)
======================

Features
--------

- A new recalculated unit cell attribute is added to the Crystal model, for use by post-integration cell refinement methods, such as that of dials.two_theta_refine. (#142)
- Add ExperimentList.change_basis() convenience method. (#166)
- Allow creation of Format classes that accept URLs instead of files (#173)


Bugfixes
--------

- Fix a bug whereby reading a single-image data set from an Eiger detector would lead to an error. (#156)
- Fix formatting of unit cell parameters with negligible standard uncertainties (#165)
- New Eiger FileWriter (20.1.16.56035) produces NeXus compliant files, which exposed a bug in finding axis sample depends on, now fixed. (#168)


Misc
----

- #164


DIALS 2.2 (2020-03-15)
======================

Bugfixes
--------

- Fix spot-finding on images with file names ending in '0000.cbf' (#133)
- Fixed imageset slicing for image sets starting from image 0 (#141)


DIALS 2.1 (2019-12-16)
======================

Features
--------

- With changes in dials.import sequences of stills are imported as individual 
  experiments all dereferencing one image set - this is the change set to support
  this on load. (#118)


Bugfixes
--------

- Reinstate support for historic VMXi EIGER 1 images (#119)
- Fix crash when opening dataset containing many .h5 files (#126)


Deprecations and Removals
-------------------------

- dxtbx extensions can no longer be imported from `dxtbx`
  and must now be imported from `dxtbx.ext` (#29)


Misc
----

- #124


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
