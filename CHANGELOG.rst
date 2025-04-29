dxtbx 3.23.0 (2025-01-08)
=========================

Features
--------

- Nexus support: Handle reading new scale_factor fields (used for detector gain). (`#756 <https://github.com/cctbx/dxtbx/issues/756>`_)
- ``dials.import``: Add a progress bar, so that it doesn't look like progress has stopped with large collections of images. (`#768 <https://github.com/cctbx/dxtbx/issues/768>`_)
- Add ``FormatSMVADSCCetaD`` to allow easier processing of 3DED images from the Ceta-D detector, which have been converted to SMV. (`#770 <https://github.com/cctbx/dxtbx/issues/770>`_)


Bugfixes
--------

- ``dials.show``: Hide progress bar if DIALS_NOBANNER (`#774 <https://github.com/cctbx/dxtbx/issues/774>`_)


Deprecations and Removals
-------------------------

- Python 3.10 is now the minimum required (`#769 <https://github.com/cctbx/dxtbx/issues/769>`_)


Misc
----

- `#767 <https://github.com/cctbx/dxtbx/issues/767>`_, `#773 <https://github.com/cctbx/dxtbx/issues/773>`_, `#775 <https://github.com/cctbx/dxtbx/issues/775>`_


dxtbx 3.22.0 (2024-10-15)
=========================

Features
--------

- Add format class to read data from the NMX ESS detector. (`#764 <https://github.com/cctbx/dxtbx/issues/764>`_)


Bugfixes
--------

- ``dxtbx.dlsnxs2cbf``: Add work around for issues with data recorded at 32-bit. (`#759 <https://github.com/cctbx/dxtbx/issues/759>`_)
- Auxiliary data processing files (mask, gain, pedestal, and dx and dy maps) will now always be loaded when available. (`#760 <https://github.com/cctbx/dxtbx/issues/760>`_)
- Allow triangles in polygon masking. (`#761 <https://github.com/cctbx/dxtbx/issues/761>`_)
- Refactor panel positions of FormatISISSXD to account for differences in panel positions, depending on the date of data collection. (`#762 <https://github.com/cctbx/dxtbx/issues/762>`_)
- Raise a more suitable error message when failing to load an experiment list. (`#763 <https://github.com/cctbx/dxtbx/issues/763>`_)


Misc
----

- `#753 <https://github.com/cctbx/dxtbx/issues/753>`_, `#754 <https://github.com/cctbx/dxtbx/issues/754>`_, `#755 <https://github.com/cctbx/dxtbx/issues/755>`_, `#758 <https://github.com/cctbx/dxtbx/issues/758>`_


Dxtbx 3.22.0 (2024-10-15)
=========================

Features
--------

- Add format class to read data from the NMX ESS detector. (`#764 <https://github.com/cctbx/dxtbx/issues/764>`_)


Bugfixes
--------

- ``dxtbx.dlsnxs2cbf``: add bit_depth_image explicitly to work around issues with data recorded at 32 bit (`#759 <https://github.com/cctbx/dxtbx/issues/759>`_)
- Ensure that data processing auxililary files (mask, gain, pedestal, and
  dx and dy maps) are loaded whenever available. This fixes
  https://github.com/dials/dials/issues/2744 (`#760 <https://github.com/cctbx/dxtbx/issues/760>`_)
- + allow triangles in polygon masking (`#761 <https://github.com/cctbx/dxtbx/issues/761>`_)
- Refactor panel positions of FormatISISSXD to account for differences in panel positions depending on the date of data collection. (`#762 <https://github.com/cctbx/dxtbx/issues/762>`_)
- Raise a more suitable error message when failing to load an experiment list. (`#763 <https://github.com/cctbx/dxtbx/issues/763>`_)


Misc
----

- `#753 <https://github.com/cctbx/dxtbx/issues/753>`_, `#754 <https://github.com/cctbx/dxtbx/issues/754>`_, `#755 <https://github.com/cctbx/dxtbx/issues/755>`_, `#758 <https://github.com/cctbx/dxtbx/issues/758>`_


DIALS 3.21.1 (2024-08-23)
=========================

Bugfixes
--------

- Stop ``dxtbx.image_average`` shuffling panel positions for segmented detectors. (`#752 <https://github.com/cctbx/dxtbx/issues/752>`_)


dxtbx 3.21.0 (2024-08-20)
=========================

Features
--------

- Add Nonius KappaCCD format. (`#741 <https://github.com/cctbx/dxtbx/issues/741>`_)


Bugfixes
--------

- ``FormatMRC``: Relax restrictive check on the overloaded MZ header value, which caused failures to read files where MZ == 1. (`#740 <https://github.com/cctbx/dxtbx/issues/740>`_)
- ``FormatCBFMini``: When parsing header lines for a timestamp, avoid clashes with Windows paths. (`#742 <https://github.com/cctbx/dxtbx/issues/742>`_)
- ``FormatPy``: Add fix for pickle files. (`#744 <https://github.com/cctbx/dxtbx/issues/744>`_)
- ``FormatSMVRigakuSaturnNoTS``: Fix a bug in reading the image pedestal from headers. (`#746 <https://github.com/cctbx/dxtbx/issues/746>`_)


Misc
----

- `#739 <https://github.com/cctbx/dxtbx/issues/739>`_, `#743 <https://github.com/cctbx/dxtbx/issues/743>`_, `#748 <https://github.com/cctbx/dxtbx/issues/748>`_, `#749 <https://github.com/cctbx/dxtbx/issues/749>`_, `#750 <https://github.com/cctbx/dxtbx/issues/750>`_


dxtbx 3.20.0 (2024-06-19)
=========================

Features
--------

- The template handling mechanism is extended so that a template with a
  single ``#`` is expanded to match non-zero padded sequential numbers.
  For example, ``image_#.cbf`` will match ``image_1.cbf``, ``image_2.cbf``,
  ..., ``image_10.cbf`` and so on.

  Using a single ``#`` to match up to 10 images _within_ a zero-padded
  sequence continues to work as before. For example,
  ``dials.import template=insulin_1_01#.img`` will match the files
  ``insulin_1_010.img``, ``insulin_1_011.img``, ..., ``insulin_1_019.img``,
  and no others. (`#705 <https://github.com/cctbx/dxtbx/issues/705>`_)
- Allows stepping through XTC streams at specific indices provided by a text file. (`#709 <https://github.com/cctbx/dxtbx/issues/709>`_)
- Compatibility with Python 3.12. (`#725 <https://github.com/cctbx/dxtbx/issues/725>`_)
- Add ``dxtbx.any2nexus`` program, to convert any file dxtbx can read to a NeXus file. (`#735 <https://github.com/cctbx/dxtbx/issues/735>`_)


Bugfixes
--------

- ``FormatROD``: set the beam probe to "electron" for 3D ED experiments. (`#728 <https://github.com/cctbx/dxtbx/issues/728>`_)
- Raise an error if ``geometry.goniometer.axis=`` is set with a multi-axis goniometer. In that case ``geometry.goniometer.axes=`` must be set instead. (`#730 <https://github.com/cctbx/dxtbx/issues/730>`_)
- Update goniometer for `FormatISISSXD` to allow for different ways the goniometer angle can be stored. (`#731 <https://github.com/cctbx/dxtbx/issues/731>`_)
- Fix `Scan.get_property` key type. (`#734 <https://github.com/cctbx/dxtbx/issues/734>`_)


Misc
----

- `#702 <https://github.com/cctbx/dxtbx/issues/702>`_, `#721 <https://github.com/cctbx/dxtbx/issues/721>`_, `#724 <https://github.com/cctbx/dxtbx/issues/724>`_, `#726 <https://github.com/cctbx/dxtbx/issues/726>`_, `#727 <https://github.com/cctbx/dxtbx/issues/727>`_, `#732 <https://github.com/cctbx/dxtbx/issues/732>`_, `#733 <https://github.com/cctbx/dxtbx/issues/733>`_, `#738 <https://github.com/cctbx/dxtbx/issues/738>`_


DIALS 3.19.1 (2024-05-23)
=========================

Bugfixes
--------

- Fix case where old I03 Eiger nexus data (pre-2020) would fail to process. (`#737 <https://github.com/cctbx/dxtbx/issues/737>`_)


dxtbx 3.19.0 (2024-04-17)
=========================

Features
--------

- Add format reader for Jungfrau4M serial images from beamline ID29 at ESRF. (`#659 <https://github.com/cctbx/dxtbx/issues/659>`_)
- Better handle spectra calibration for bad data in XTC format using new parameter: ``spectrum_required=``. (`#674 <https://github.com/cctbx/dxtbx/issues/674>`_)
- Add Bruker and miniCBF format readers for the ELDICO ED-1 electron diffractometer with DECTRIS QUADRO detector. (`#682 <https://github.com/cctbx/dxtbx/issues/682>`_)
- ``FormatSMVTimePix_SU``: Always mask out the central cross of virtual pixels. (`#683 <https://github.com/cctbx/dxtbx/issues/683>`_)
- Add format reader for ISIS SXD detector. (`#687 <https://github.com/cctbx/dxtbx/issues/687>`_)
- Detector distance can now be manually overridden for multi-panel detectors. (`#698 <https://github.com/cctbx/dxtbx/issues/698>`_)
- Add format reader to read time of flight Laue data from MANDI. (`#703 <https://github.com/cctbx/dxtbx/issues/703>`_)
- Additional features for `FormatXTCRayonix` (`#723 <https://github.com/cctbx/dxtbx/issues/723>`_)


Bugfixes
--------

- Importing the (deprecated and removed) ``dxtbx.datablock`` module failed to display warning properly. (`#665 <https://github.com/cctbx/dxtbx/issues/665>`_)
- Fix scan comparison for scan properties changes (`#669 <https://github.com/cctbx/dxtbx/issues/669>`_)
- Eiger Support: Invert the module dimensions, only for older firmware versions. See https://media.dectris.com/230203-Release_Notes-DECTRIS_EIGER2.pdf for reference. (`#676 <https://github.com/cctbx/dxtbx/issues/676>`_)
- ``FormatMRC``: Better handling of extended headers. (https://github.com/ccpem/mrcfile/issues/50), and extended headers are ignored if they contain junk values. (`#679 <https://github.com/cctbx/dxtbx/issues/679>`_)
- Fixed some properties not being correctly parsed in `Scan.from_dict`. (`#688 <https://github.com/cctbx/dxtbx/issues/688>`_)
- Negative rotation angles are now allowed, the goniometer axis will be inverted if necessary. (`#690 <https://github.com/cctbx/dxtbx/issues/690>`_)
- ``dials.import`` now uses natural sorting on input data, instead of strict sorting. (`#697 <https://github.com/cctbx/dxtbx/issues/697>`_)
- Fix setting detector distance for single panel detectors that have a hierarchy. (`#699 <https://github.com/cctbx/dxtbx/issues/699>`_)
- Better recognition for SMV images from MLFSOM and other simulators from James Holton. (`#708 <https://github.com/cctbx/dxtbx/issues/708>`_)
- Fix error introduced in ``FormatSMVJHSim`` causing test failures. (`#710 <https://github.com/cctbx/dxtbx/issues/710>`_)
- `PolychromaticBeam` can now be copied with `copy.deepcopy`. (`#711 <https://github.com/cctbx/dxtbx/issues/711>`_)
- Add missing argument to `PolychromaticBeamPickleSuite.getinitargs`. (`#714 <https://github.com/cctbx/dxtbx/issues/714>`_)


Misc
----

- `#620 <https://github.com/cctbx/dxtbx/issues/620>`_, `#667 <https://github.com/cctbx/dxtbx/issues/667>`_, `#670 <https://github.com/cctbx/dxtbx/issues/670>`_, `#689 <https://github.com/cctbx/dxtbx/issues/689>`_, `#691 <https://github.com/cctbx/dxtbx/issues/691>`_, `#694 <https://github.com/cctbx/dxtbx/issues/694>`_, `#696 <https://github.com/cctbx/dxtbx/issues/696>`_, `#701 <https://github.com/cctbx/dxtbx/issues/701>`_, `#704 <https://github.com/cctbx/dxtbx/issues/704>`_, `#707 <https://github.com/cctbx/dxtbx/issues/707>`_, `#713 <https://github.com/cctbx/dxtbx/issues/713>`_


dxtbx 3.17.0 (2023-11-03)
=========================

Features
--------

- Add ``nxmx_writer``, a tool for converting any data dxtbx can read to NeXus data. (`#615 <https://github.com/cctbx/dxtbx/issues/615>`_)
- Remove circular dependencies between dxtbx and ``cctbx.xfel``, by using the new ``serialtbx``. (`#627 <https://github.com/cctbx/dxtbx/issues/627>`_)
- Set the beam probe to ``electron`` in both ``FormatNXmxED`` and ``FormatSER``. (`#661 <https://github.com/cctbx/dxtbx/issues/661>`_)


Bugfixes
--------

- ``dxtbx.image_average``: Better handle detector gain and pixel data type. (`#660 <https://github.com/cctbx/dxtbx/issues/660>`_)
- ``Beam.probe`` is no longer reset if any geometrical override is provided at import. (`#661 <https://github.com/cctbx/dxtbx/issues/661>`_)
- Pilatus 4: Do not invert module size that is correctly written in master file. (`#663 <https://github.com/cctbx/dxtbx/issues/663>`_)
- ``dxtbx.plot_detector_models``: Use noninteractive matpotlib backend, if using the ``pdf_file=`` option. (`#664 <https://github.com/cctbx/dxtbx/issues/664>`_)


Deprecations and Removals
-------------------------

- Legacy ``Datablock`` support has been removed, after being deprecated for several years. If you have any experiments that use these, they will need to be re-imported. (`#504 <https://github.com/cctbx/dxtbx/issues/504>`_)


Misc
----

- `#622 <https://github.com/cctbx/dxtbx/issues/622>`_


Dxtbx 3.17 (2023-11-03)
=======================

Features
--------

- Add nxmx_writer, a tool for converting any data dxtbx can read to NeXus data (`#615 <https://github.com/cctbx/dxtbx/issues/615>`_)
- Remove circular dependencies between dxtbx and ``cctbx.xfel`` by using the new ``serialtbx``. (`#627 <https://github.com/cctbx/dxtbx/issues/627>`_)
- Set the beam probe to ``electron`` in both ``FormatNXmxED`` and ``FormatSER``. (`#661 <https://github.com/cctbx/dxtbx/issues/661>`_)


Bugfixes
--------

- Bugfix for dxtbx.image_average: handle detector gain and pixel data type better (`#660 <https://github.com/cctbx/dxtbx/issues/660>`_)
- The beam probe is no longer reset if any geometrical override is provided at import. (`#661 <https://github.com/cctbx/dxtbx/issues/661>`_)
- Pilatus 4: do not invert module size (is written correctly in master file) (`#663 <https://github.com/cctbx/dxtbx/issues/663>`_)
- ``dxtbx.plot_detector_models``: use noninteractive matpotlib backend if using the pdf_file option (`#664 <https://github.com/cctbx/dxtbx/issues/664>`_)


Deprecations and Removals
-------------------------

- dxtbx: remove legacy datablock object (obsolete for several years) (`#504 <https://github.com/cctbx/dxtbx/issues/504>`_)


Misc
----

- `#622 <https://github.com/cctbx/dxtbx/issues/622>`_


DIALS 3.16.1 (2023-09-05)
=========================

Bugfixes
--------

- Fix situation where a bad ``Beam.probe`` could cause undefined behaviour. (`#656 <https://github.com/cctbx/dxtbx/issues/656>`_)
- Fix performance regression loading large experiment lists containing profile/scaling models. (`#658 <https://github.com/cctbx/dxtbx/issues/658>`_)


dxtbx 3.16.0 (2023-08-14)
=========================

Features
--------

- Add new Beam class ``dxtbx.model.PolychromaticBeam``, for polychromatic/multi-wavelength/wide bandpass experiments. (`#621 <https://github.com/cctbx/dxtbx/issues/621>`_)
- Formats: Reflect move of Eiger detector from PETRA P14 to P13. (`#626 <https://github.com/cctbx/dxtbx/issues/626>`_)
- The ``model.Beam`` object now has a ``probe`` value to keep track of the type of radiation. (`#647 <https://github.com/cctbx/dxtbx/issues/647>`_)
- Formats: CBFMini support for the EIGER2 16M detector at CHESS beamline ID7B2, which has an inverted rotation axis. (`#649 <https://github.com/cctbx/dxtbx/issues/649>`_)
- Formats: Support for Eiger 9M on ESRF ID23-2, which has an undeclared vertical goniometer. (`#651 <https://github.com/cctbx/dxtbx/issues/651>`_)
- Formats: Partial support for the Rigaku Oxford Diffraction file format, including support for multi-axis goniometers and faster decompression. (`#645 <https://github.com/cctbx/dxtbx/issues/645>`_) (`#653 <https://github.com/cctbx/dxtbx/issues/653>`_)


Bugfixes
--------

- Panel geometry definitions in PHIL are merged by panel id *before* constructing panels. (`#299 <https://github.com/cctbx/dxtbx/issues/299>`_)
- ``flumpy``: Fix case where incorrect ``flex.vec2``, ``flex.vec3`` could be generated. (`#439 <https://github.com/cctbx/dxtbx/issues/439>`_)
- NXmx files with multidimensional arrays (images, modules, or both) are now handled. (`#612 <https://github.com/cctbx/dxtbx/issues/612>`_)
- Slicing of imageset objects is now consistently 0-based, including for the sliced data accessor. Previously, the data accessor had to be accessed with the original index offsets. (`#633 <https://github.com/cctbx/dxtbx/issues/633>`_)
- Formats: Add fix for Eiger / NXmx data from DLS i19-2, to correctly assign the image bit depth. (`#652 <https://github.com/cctbx/dxtbx/issues/652>`_)


Misc
----

- `#640 <https://github.com/cctbx/dxtbx/issues/640>`_, `#642 <https://github.com/cctbx/dxtbx/issues/642>`_, `#643 <https://github.com/cctbx/dxtbx/issues/643>`_, `#645 <https://github.com/cctbx/dxtbx/issues/645>`_, `#650 <https://github.com/cctbx/dxtbx/issues/650>`_, `#655 <https://github.com/cctbx/dxtbx/issues/655>`_


DIALS 3.15.1 (2023-06-29)
=========================

Bugfixes
--------

- ``dxtbx.dlsnxs2cbf``: Fix import overwritten by local variable. (`#641 <https://github.com/cctbx/dxtbx/issues/641>`_)


dxtbx 3.15.0 (2023-06-13)
=========================

Features
--------

- Support for Bruker Photon detectors has been extended to include Photon-III. (`#637 <https://github.com/cctbx/dxtbx/issues/637>`_)


Bugfixes
--------

- Rigaku Saturn SMV images with multi-axis crystal goniometers are now handledi, instead of being silently ignored. With thanks to James Hester for this contribution. (`#617 <https://github.com/cctbx/dxtbx/issues/617>`_)
- FormatCBFFull: If rotation angles are decreasing, then invert the rotation axis as well as the angles, to be consistent. (`#623 <https://github.com/cctbx/dxtbx/issues/623>`_)
- Bugfix for CCTBX bootstrapped environments, without conda. (`#630 <https://github.com/cctbx/dxtbx/issues/630>`_)


Misc
----

- `#625 <https://github.com/cctbx/dxtbx/issues/625>`_, `#636 <https://github.com/cctbx/dxtbx/issues/636>`_, `#639 <https://github.com/cctbx/dxtbx/issues/639>`_


DIALS 3.14.2 (2023-05-16)
=========================

Bugfixes
--------

- Compatibility fix for the DECTRIS Eiger FileWriter. Recent FileWriter versions split bit depth metadata into two separate items, ``bit_depth_readout`` from the NXmx standard, and the new ``bit_depth_image`` field. This adds support for the latter, and now passes the metadata through into image conversion. (`#632 <https://github.com/cctbx/dxtbx/issues/632>`_)


dxtbx 3.14.0 (2023-04-12)
=========================

Features
--------

- ``flumpy``: Add support for conversion of ``flex.miller_index`` arrays to/from numpy. (`#618 <https://github.com/cctbx/dxtbx/issues/618>`_)


Bugfixes
--------

- Flumpy: Prefer returning ``flex.int`` instead of ``flex.long`` when they are the same size. This solves ambiguous behaviour when reading images on Windows platforms. (`#607 <https://github.com/cctbx/dxtbx/issues/607>`_)
- ``dxtbx.plot_detector_models``: Fix display of multiple single-panel detector models. (`#610 <https://github.com/cctbx/dxtbx/issues/610>`_)


Misc
----

- `#604 <https://github.com/cctbx/dxtbx/issues/604>`_, `#608 <https://github.com/cctbx/dxtbx/issues/608>`_, `#609 <https://github.com/cctbx/dxtbx/issues/609>`_, `#611 <https://github.com/cctbx/dxtbx/issues/611>`_, `#614 <https://github.com/cctbx/dxtbx/issues/614>`_


dxtbx 3.13.0 (2023-01-26)
=========================

Features
--------

- ``FormatNXmxED``: Format support for electron diffraction images converted to be compatible with NXmx by `nexgen <https://github.com/dials/nexgen>`_. (`#583 <https://github.com/cctbx/dxtbx/issues/583>`_)
- ``FormatNXmxEDeBIC``: Including a mask specific for the SINGLA that is temporarily installed at eBIC, through to early 2023. (`#589 <https://github.com/cctbx/dxtbx/issues/589>`_)


Bugfixes
--------

- ``dxtbx.image_average``: Fix a crash from using more processors than images when using MPI. (`#571 <https://github.com/cctbx/dxtbx/issues/571>`_)
- ``dxtbx.plot_detector_models`` now works with newer versions of matplotlib. (`#574 <https://github.com/cctbx/dxtbx/issues/574>`_)
- ``FormatNXmxDLS``: Don't process electron diffraction images collected at eBIC that have been converted by ``nexgen``. (`#579 <https://github.com/cctbx/dxtbx/issues/579>`_)
- Correct maximum value of Rayonix trusted range. (`#590 <https://github.com/cctbx/dxtbx/issues/590>`_)
- Read underload from CBF files. (`#592 <https://github.com/cctbx/dxtbx/issues/592>`_)
- ``NXmx``: Ensure integer data types get converted to ``flex.int`` on all platforms. (`#594 <https://github.com/cctbx/dxtbx/issues/594>`_)
- Fix trusted range in ``FormatCBFMultiTile`` and ``FormatCBFMultiTileHierarchy``. (`#595 <https://github.com/cctbx/dxtbx/issues/595>`_)
- ``FullCBFWriter``: Fix writing of the newly consistent trusted_range values. (`#601 <https://github.com/cctbx/dxtbx/issues/601>`_)


Misc
----

- `#578 <https://github.com/cctbx/dxtbx/issues/578>`_, `#591 <https://github.com/cctbx/dxtbx/issues/591>`_, `#597 <https://github.com/cctbx/dxtbx/issues/597>`_, `#598 <https://github.com/cctbx/dxtbx/issues/598>`_, `#599 <https://github.com/cctbx/dxtbx/issues/599>`_, `#600 <https://github.com/cctbx/dxtbx/issues/600>`_, `#602 <https://github.com/cctbx/dxtbx/issues/602>`_, `#603 <https://github.com/cctbx/dxtbx/issues/603>`_, `#605 <https://github.com/cctbx/dxtbx/issues/605>`_, `#606 <https://github.com/cctbx/dxtbx/issues/606>`_


Dxtbx 3.13 (2023-01-12)
=======================

Features
--------

- ``FormatNXmxED``: Format support for electron diffraction images converted to be compatible with NXmx by `nexgen <https://github.com/dials/nexgen>`_. (`#583 <https://github.com/cctbx/dxtbx/issues/583>`_)
- ``FormatNXmxEDeBIC``: Including a mask specific for the SINGLA that is temporarily installed at eBIC, through to early 2023. (`#589 <https://github.com/cctbx/dxtbx/issues/589>`_)


Bugfixes
--------

- ``dxtbx.image_average``: Fix a crash from using more processors than images when using MPI. (`#571 <https://github.com/cctbx/dxtbx/issues/571>`_)
- ``dxtbx.plot_detector_models`` now works with newer versions of matplotlib. (`#574 <https://github.com/cctbx/dxtbx/issues/574>`_)
- ``FormatNXmxDLS``: Don't process electron diffraction images collected at eBIC that have been converted by ``nexgen``. (`#579 <https://github.com/cctbx/dxtbx/issues/579>`_)
- Correct maximum value of Rayonix trusted range. (`#590 <https://github.com/cctbx/dxtbx/issues/590>`_)
- Read underload from CBF files (`#592 <https://github.com/cctbx/dxtbx/issues/592>`_)
- ``NXmx``: Ensure integer data types get converted to ``flex.int`` on all platforms. (`#594 <https://github.com/cctbx/dxtbx/issues/594>`_)
- Fix trusted range in ``FormatCBFMultiTile`` and ``FormatCBFMultiTileHierarchy``. (`#595 <https://github.com/cctbx/dxtbx/issues/595>`_)


Misc
----

- `#578 <https://github.com/cctbx/dxtbx/issues/578>`_, `#591 <https://github.com/cctbx/dxtbx/issues/591>`_, `#597 <https://github.com/cctbx/dxtbx/issues/597>`_, `#598 <https://github.com/cctbx/dxtbx/issues/598>`_, `#600 <https://github.com/cctbx/dxtbx/issues/600>`_


Dxtbx 3.13 (2023-01-12)
=======================

Features
--------

- ``FormatNXmxED``: a new format class for electron diffraction images converted to be compatible with NXmx by nexgen (https://github.com/dials/nexgen) (`#583 <https://github.com/cctbx/dxtbx/issues/583>`_)
- Add ``FormatNXmxEDeBIC``, which includes a mask specific for the SINGLA that is temporarily installed at eBIC, through to early 2023. (`#589 <https://github.com/cctbx/dxtbx/issues/589>`_)


Bugfixes
--------

- dxtbx.image_average: fix a crash from using more processors than images when using MPI. (`#571 <https://github.com/cctbx/dxtbx/issues/571>`_)
- Fix dxtbx.plot_detector_models for new versions of matplotlib (`#574 <https://github.com/cctbx/dxtbx/issues/574>`_)
- ``FormatNXmxDLS`` no longer recognises electron diffraction images collected at eBIC that have been converted by ``nexgen``. (`#579 <https://github.com/cctbx/dxtbx/issues/579>`_)
- Corrected maximum value of Rayonix trusted range. (`#590 <https://github.com/cctbx/dxtbx/issues/590>`_)
- Read underload from CBF files (`#592 <https://github.com/cctbx/dxtbx/issues/592>`_)
- ``NXmx``: Ensure integer data types get converted to ``flex.int``, i.e. ``int`` C-type, on all platforms (`#594 <https://github.com/cctbx/dxtbx/issues/594>`_)
- Fix trusted range in FormatCBFMultiTile and FormatCBFMultiTileHierarchy. (`#595 <https://github.com/cctbx/dxtbx/issues/595>`_)


Misc
----

- `#578 <https://github.com/cctbx/dxtbx/issues/578>`_, `#591 <https://github.com/cctbx/dxtbx/issues/591>`_, `#597 <https://github.com/cctbx/dxtbx/issues/597>`_, `#598 <https://github.com/cctbx/dxtbx/issues/598>`_, `#600 <https://github.com/cctbx/dxtbx/issues/600>`_


DIALS 3.12.1 (2022-12-05)
=========================

Bugfixes
--------

- NXmx: eliminate potential divide-by-zero warning (`#572 <https://github.com/cctbx/dxtbx/issues/572>`_)
- Fallback on legacy FormatNexus to workaround issues reading datasets written by the Dectris filewriter with FormatNXmx(#582) (`#584 <https://github.com/cctbx/dxtbx/issues/584>`_)
- Fix support for datasets generated by the DECTRIS EIGER filewriter (`#586 <https://github.com/cctbx/dxtbx/issues/586>`_)
- ``FormatCBFFull``: trusted range bug fix - use the minimum valid pixel value rather than the undefined value (`#587 <https://github.com/cctbx/dxtbx/issues/587>`_)
- NXmx: fallback on explicit beam_center_{x,y} if the x,y components of the detector origin are zero (`#588 <https://github.com/cctbx/dxtbx/issues/588>`_)


dxtbx 3.12.0 (2022-10-31)
=========================

Features
--------

- Improve XTC handling from LCLS. Includes better spectrum support, parallax for the ePix, binning for the Rayonix, and radial_averge fixes. (`#517 <https://github.com/cctbx/dxtbx/issues/517>`_)
- Add spectrum support to FormatNXmx. (`#538 <https://github.com/cctbx/dxtbx/issues/538>`_)
- NXmx: Add support for `@equipment_component <https://manual.nexusformat.org/classes/base_classes/NXtransformations.html#nxtransformations-axisname-equipment-component-attribute>`_ for forming logical groupings of transformations to reduce the number of levels in the detector hierarchy.  Note: ``.expt`` files will not be backwards compatible for users of the JF16M detector at SwissFEL, or the AGIPD detector at EuXFEL. (`#561 <https://github.com/cctbx/dxtbx/issues/561>`_)


Bugfixes
--------

- ``trusted_range`` is now defined consistently as the _inclusive_ range between the minimum and maximum trusted values, i.e. valid pixels are those less than or equal to the maximum trusted value and greater than or equal to the minimum trusted value. (`#536 <https://github.com/cctbx/dxtbx/issues/536>`_)
- Improved speed of reading many-panel Nexus images. (`#565 <https://github.com/cctbx/dxtbx/issues/565>`_)
- Remove unintended error message escalation when passing multiple image ranges to import. 
- Remove stray and unhelpful error message display when passing multiple image ranges to import. (`#567 <https://github.com/cctbx/dxtbx/issues/567>`_)
- Added Diamonds VMXm Eiger CdTe 9M to "legacy" list where the fast, slow dimensions are reversed. (`#569 <https://github.com/cctbx/dxtbx/issues/569>`_)


Deprecations and Removals
-------------------------

- The deprecated ``set_slow_fast_beam_centre_mm`` function has been removed. Please use ``set_fast_slow_beam_centre_mm`` instead. (`#544 <https://github.com/cctbx/dxtbx/issues/544>`_)


Misc
----

- `#541 <https://github.com/cctbx/dxtbx/issues/541>`_, `#543 <https://github.com/cctbx/dxtbx/issues/543>`_, `#554 <https://github.com/cctbx/dxtbx/issues/554>`_, `#556 <https://github.com/cctbx/dxtbx/issues/556>`_, `#557 <https://github.com/cctbx/dxtbx/issues/557>`_, `#558 <https://github.com/cctbx/dxtbx/issues/558>`_, `#563 <https://github.com/cctbx/dxtbx/issues/563>`_


DIALS 3.11.2 (2022-09-27)
=========================

Bugfixes
--------

- ``NXmx``: Cope more gracefully with scalar NXtransformations values. (`#546 <https://github.com/cctbx/dxtbx/issues/546>`_)
- ``dxtbx.dlsnxs2cbf``: Fix distance and pixel size bugs. (`#548 <https://github.com/cctbx/dxtbx/issues/548>`_)
- NXmx reading: Handle cases where the detector is read as between the sample and source. This is to compensate for an incorrect definition in the Dectris Eiger file writer. (`#550 <https://github.com/cctbx/dxtbx/issues/550>`_)


Misc
----

- `#547 <https://github.com/cctbx/dxtbx/issues/547>`_


DIALS 3.11.1 (2022-09-02)
=========================

Bugfixes
--------

- ``dxtbx.dlsnxs2cbf``: Fix bug introduced by #572. (`#545 <https://github.com/cctbx/dxtbx/issues/545>`_)


dxtbx 3.11.0 (2022-08-24)
=========================

Features
--------

- Replace use of legacy ``FormatNexusEiger`` with new ``FormatNXmx`` format class. (`#455 <https://github.com/cctbx/dxtbx/issues/455>`_)


Bugfixes
--------

- DXTBX now uses the median oscillation width from across the entire scan. This resolved issues where the goniometer scan positions were read-back values instead of set-point values, and a slow rotation start across the first two images would cause the oscillation width for the whole scan to be calculated incorrectly. (`#526 <https://github.com/cctbx/dxtbx/issues/526>`_)
- ``FormatNXmx``: Support NXmx files with one wavelength per image. (`#527 <https://github.com/cctbx/dxtbx/issues/527>`_)
- ``ExperimentList.append()``: No longer O(N²) with experiment identifiers. (`#528 <https://github.com/cctbx/dxtbx/issues/528>`_)
- ``FormatNXmx``: Ignore empty pixel masks, instead of printing a confusing error. (`#529 <https://github.com/cctbx/dxtbx/issues/529>`_)
- Correct assumptions about interpreting multi-axis goniometer axes from full-CBF files. Previously, it was assumed the ``axis`` and ``diffrn_scan_axis`` categories listed axes in the same order, and that this matched a standard diffractometer axis order. The goniometer model is now build correctly, regardless of the order specified in the file. (`#539 <https://github.com/cctbx/dxtbx/issues/539>`_)


Misc
----

- `#531 <https://github.com/cctbx/dxtbx/issues/531>`_, `#533 <https://github.com/cctbx/dxtbx/issues/533>`_


DIALS 3.10.3 (2022-08-02)
=========================

Bugfixes
--------

- Fix ``mask_untrusted_circle()`` crash when untrusted circle extends outside detector. This affected ``dials.generate_mask``. (`#525 <https://github.com/cctbx/dxtbx/issues/525>`_)
- ``FormatNXmx``: Allow empty ``saturation_value`` field when importing data. (`#534 <https://github.com/cctbx/dxtbx/issues/534>`_)


DIALS 3.10.1 (2022-07-12)
=========================

Features
--------

- Updated bad pixel mask for DLS I23 PILATUS 12M for 2022 run 3 (`#530 <https://github.com/cctbx/dxtbx/issues/530>`_)


Bugfixes
--------

- ``dxtbx.install_format``: Handle case on MacOS ``.pkg`` installations where URL-formats could not be installed. (`#524 <https://github.com/cctbx/dxtbx/issues/524>`_)


dxtbx 3.10.0 (2022-06-09)
=========================

Features
--------

- Recognise `NXmx standard <https://manual.nexusformat.org/classes/applications/NXmx.html>`_ data from the Diamond Light Source `DIAD <https://www.diamond.ac.uk/Instruments/Imaging-and-Microscopy/DIAD.html>`_ beamline. (`#506 <https://github.com/cctbx/dxtbx/issues/506>`_)
- When installed as a libtbx module, dxtbx will not install python packages into ``conda_base/``. (`#511 <https://github.com/cctbx/dxtbx/issues/511>`_)
- Added ``flex_table.h`` and ``flex_table_suite.h`` objects from DIALS. These contain the C++ classes backing the ``dials.array_family.flex.reflection_table`` object, and allow a collection of ``array_family.flex`` arrays to be grouped together into a multi-columnar, row-addressable format. They are moved here to allow extension of the dxtbx models in this form. (`#521 <https://github.com/cctbx/dxtbx/issues/521>`_)


Bugfixes
--------

- Fixed ``Panel.projection_2d`` not being serialized. (`#509 <https://github.com/cctbx/dxtbx/issues/509>`_)
- ``dxtbx.dlsnxs2cbf``: Fix image oscillation for screening images (`#514 <https://github.com/cctbx/dxtbx/issues/514>`_)
- Fix ``dxtbx.image_average`` for raster scans. (`#522 <https://github.com/cctbx/dxtbx/issues/522>`_)


Deprecations and Removals
-------------------------

- Remove disused ``FormatEigerStream`` format class. This was used internally at Diamond Light Source as an intermediate solution before implementing SWMR support. (`#499 <https://github.com/cctbx/dxtbx/issues/499>`_)


Misc
----

- `#498 <https://github.com/cctbx/dxtbx/issues/498>`_, `#500 <https://github.com/cctbx/dxtbx/issues/500>`_, `#502 <https://github.com/cctbx/dxtbx/issues/502>`_, `#505 <https://github.com/cctbx/dxtbx/issues/505>`_, `#512 <https://github.com/cctbx/dxtbx/issues/512>`_, `#513 <https://github.com/cctbx/dxtbx/issues/513>`_, `#515 <https://github.com/cctbx/dxtbx/issues/515>`_, `#520 <https://github.com/cctbx/dxtbx/issues/520>`_


dxtbx DIALS 3.9.2 (2022-05-09)
==============================

Bugfixes
--------

- ``FormatCBFFullPilatus``: Handle detector information better in cases of multiple or missing panels. (`#508 <https://github.com/cctbx/dxtbx/issues/508>`_)
- Remove check for beam/normalization orthogonality in ``Beam.rotate_around_origin``. This could stop processing of older, incorrectly configured data. (`#510 <https://github.com/cctbx/dxtbx/issues/510>`_)
- Correct a unicode error reading Bruker ``.sfrm`` files. With thanks to `Dennis Brookner <https://github.com/dennisbrookner>`_ for this change. (`#518 <https://github.com/cctbx/dxtbx/issues/518>`_)


dxtbx 3.8.4 (2022-04-01)
========================

Bugfixes
--------

- ``FormatNXmxI19_2``:  Allow data from beamline I19-2 at Diamond Light Source to be processed with optional masking of the beamline's standard diamond anvil pressure cell with a 76° aperture. (`#481 <https://github.com/cctbx/dxtbx/issues/481>`_)


dxtbx 3.9.1 (2022-03-31)
========================

Features
--------

- Windows support for the CMake build. (`#507 <https://github.com/cctbx/dxtbx/issues/507>`_)


dxtbx 3.9.0 (2022-03-14)
========================

Features
--------

- Add get_spectrum to FormatXTC (`#484 <https://github.com/cctbx/dxtbx/issues/484>`_)
- Add filtering by event code for processing LCLS data (`#489 <https://github.com/cctbx/dxtbx/issues/489>`_)
- Beam flux is now written to, and read from, CBF files. (`#493 <https://github.com/cctbx/dxtbx/issues/493>`_)


Bugfixes
--------

- Reduce, in some cases drastically, memory usage of ``ImageSet`` objects. (`#438 <https://github.com/cctbx/dxtbx/issues/438>`_)
- Make FormatPY abstract so that dxtbx doesn't try to read ``.pickle`` reflection files as images. (`#464 <https://github.com/cctbx/dxtbx/issues/464>`_)
- Add method ersatz_uuid4 which gives an implementation of a random 128 bit UUID4 (`#477 <https://github.com/cctbx/dxtbx/issues/477>`_)
- ``FormatNXmxI19_2``:  Allow data from beamline I19-2 at Diamond Light Source to be processed with optional masking of the beamline's standard diamond anvil pressure cell with a 76° aperture. (`#481 <https://github.com/cctbx/dxtbx/issues/481>`_)
- Correctly handle slicing ImageSequences made from images starting with 0 (`#485 <https://github.com/cctbx/dxtbx/issues/485>`_)
- The Beam object constructor no longer discards "transmission" and "flux". (`#488 <https://github.com/cctbx/dxtbx/issues/488>`_)
- Fix wavelength bug in FormatXTC for older datasets (`#490 <https://github.com/cctbx/dxtbx/issues/490>`_)
- Fixed inconsistency in ``dxtbx.model.Scan`` default constructor that gave different results when loading from Python dictionary. (`#496 <https://github.com/cctbx/dxtbx/issues/496>`_)


Misc
----

- `#462 <https://github.com/cctbx/dxtbx/issues/462>`_, `#463 <https://github.com/cctbx/dxtbx/issues/463>`_, `#466 <https://github.com/cctbx/dxtbx/issues/466>`_, `#468 <https://github.com/cctbx/dxtbx/issues/468>`_, `#471 <https://github.com/cctbx/dxtbx/issues/471>`_, `#477 <https://github.com/cctbx/dxtbx/issues/477>`_, `#479 <https://github.com/cctbx/dxtbx/issues/479>`_, `#480 <https://github.com/cctbx/dxtbx/issues/480>`_, `#482 <https://github.com/cctbx/dxtbx/issues/482>`_, `#487 <https://github.com/cctbx/dxtbx/issues/487>`_, `#494 <https://github.com/cctbx/dxtbx/issues/494>`_, `#495 <https://github.com/cctbx/dxtbx/issues/495>`_


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

- dxtbx can be optionally used without ``cbflib_adaptbx``. (`#368 <https://github.com/cctbx/dxtbx/issues/368>`_)
- Experimental support for building dxtbx with CMake. (`#449 <https://github.com/cctbx/dxtbx/issues/449>`_)
- Track dxtbx version explicitly, with bump2version. (`#458 <https://github.com/cctbx/dxtbx/issues/458>`_)


Bugfixes
--------

- Fix an arithmetic mistake in ``dxtbx.model.Goniometer.rotate_around_origin``, which was mangling the addition of a new rotation to the goniostat rotation operator :math:`\mathbf{R}`. (`#451 <https://github.com/cctbx/dxtbx/issues/451>`_)
- Correct pedestal handling for simulated images from ``simtbx``. (`#456 <https://github.com/cctbx/dxtbx/issues/456>`_)
- Ensure ``FormatTIFF`` only understands images with the expected basic TIFF header. (`#457 <https://github.com/cctbx/dxtbx/issues/457>`_)
- Get CI builds working again by restricting ``setuptools<60``. (`#459 <https://github.com/cctbx/dxtbx/issues/459>`_)


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
