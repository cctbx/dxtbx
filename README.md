## Diffraction Experiment Toolbox

![Python 2.7 | 3.5 | 3.6 | 3.7](https://img.shields.io/badge/python-2.7%20%7C%203.5%20%7C%203.6%20%7C%203.7-blue.svg)
![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/dials/dxtbx.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/dials/dxtbx/context:python)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/dials/dxtbx.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/dials/dxtbx/alerts/)
[![Coverage](https://codecov.io/github/dials/dxtbx/coverage.svg?branch=master)](https://codecov.io/github/dials/dxtbx?branch=master)

A [cctbx](https://cctbx.github.io/)-style toolbox to describe single-crystal diffraction experiments, where
a monochromatic beam is used to illuminate a sample which is rotated during
the exposure and diffraction recorded on a flat area detector.

This toolbox will include code for:

 * reading image headers
 * transforming contents of image header to standard (i.e. imgCIF) frame
 * python models of experiment
 * reading a sequence into memory using existing cctbx image reading tools in [iotbx](https://cctbx.github.io/iotbx/index.html)

Initially implemented to support [xia2](https://github.com/xia2/xia2) development, dxtbx is designed to be extensible, to support other applications and to make it easy to work with other detectors, with a generic approach to reading the data files.

A paper describing how to use dxtbx, as well as documenting its development and some of its applications, was published as [J. Appl. Cryst. (2014) **47**, 1459-1465](https://doi.org/10.1107/S1600576714011996).
