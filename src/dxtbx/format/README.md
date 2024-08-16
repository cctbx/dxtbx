## dxtbx format classes

dxtbx uses a collection of format classes to understand image data.
To add support for a new image format a new format class is required.
All format classes inherit directly or indirectly from a common `Format`
class, and by convention, all format classes have a name beginning with 'Format'.

Format classes do not have to be distributed as part of dxtbx, they
can be distributed separately. Any python package can declare contained
format classes using entry points in setup.py as follows:

```python
setup(
    ...,
    entry_points={
        "dxtbx.format": [
            "CLASSNAME:PARENTS = python.module.path:CLASSNAME"
        ],
    },
    ...,
}
```
For example to declare a format class FormatMyClass from file
/subdir/myformat.py in python package 'mypackage', which is a subclass
of both FormatCBF and FormatMultiImage, the full entry point string is:

    FormatMyClass:FormatCBF,FormatMultiImage = mypackage.subdir.myformat:FormatMyClass

cctbx/libtbx packages do not follow the python setup.py convention, but they
can declare entry points in `libtbx_refresh.py` using
```python
libtbx.pkg_utils.define_entry_points({"dxtbx.format": [ ... ]})
```
(Note that there should only be a single `define_entry_points` call per
`libtbx_refresh.py`, as later invocations overwrite previous declarations.)

Alternatively, format classes can be distributed outside of python packages
and installed by the user using `dxtbx.install_format`. For the correct identification
of the format class and its dependencies the format class must be declared
unconditionally and the parent format class names must be imported into the
global namespace:
```python
if something:
    def FormatClass(Format):  # conditional declaration, invalid

def FormatClass(path.to.module.FormatParent):  # parent not in global namespace, invalid

from path.to.module import FormatParent1, FormatParent2
def FormatClass(FormatParent1, FormatParent2):  # valid
```
For more information on installing format classes this way please check
```
dxtbx.install_format --help
```

