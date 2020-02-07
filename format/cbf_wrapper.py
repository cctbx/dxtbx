from __future__ import print_function
import pycbf
from libtbx.utils import to_bytes, to_str


class cbf_wrapper(pycbf.cbf_handle_struct):
    """Wrapper class that provides convenience functions for working with cbflib.
       It ensures that string arguments are passed to parent class as bytes,
       while bytes returned from the parent class are converted to strings."""

    def add_category(self, name, columns):
        """Create a new category and populate it with column names"""
        self.new_category(to_bytes(name))
        for column in columns:
            self.new_column(to_bytes(column))

    def add_row(self, data):
        """Add a row to the current category.  If data contains more entries than
        there are columns in this category, then the remainder is truncated
        Use '.' for an empty value in a row."""

        self.new_row()
        self.rewind_column()
        for item in data:
            self.set_value(to_bytes(item))
            if item == ".":
                self.set_typeofvalue(b"null")
            try:
                self.next_column()
            except Exception:
                break

    def has_sections(self):
        """True if the cbf has the array_structure_list_section table, which
        changes how its data is stored in the binary sections
        """
        try:
            self.find_category("array_structure_list_section")
            return True
        except Exception as e:
            if "CBF_NOTFOUND" in str(e):
                return False
            raise e

    def read_widefile(self, imagefile, flags=pycbf.MSG_DIGEST):
        return super(cbf_wrapper, self).read_widefile(to_bytes(imagefile), flags)

    def new_datablock(self, arg):
        return to_str(super(cbf_wrapper, self).new_datablock(to_bytes(arg)))

    def find_category(self, arg):
        return to_str(super(cbf_wrapper, self).find_category(to_bytes(arg)))

    def find_column(self, arg):
        return to_str(super(cbf_wrapper, self).find_column(to_bytes(arg)))

    def set_realarray_wdims_fs(
        self,
        pycbf_const,
        binary_id,
        data,
        elsize,
        elements,
        byteorder,
        dimfast,
        dimmid,
        dimslow,
        padding,
    ):
        return super(cbf_wrapper, self).set_realarray_wdims_fs(
            pycbf_const,
            binary_id,
            data,
            elsize,
            elements,
            to_bytes(byteorder),
            dimfast,
            dimmid,
            dimslow,
            padding,
        )

    def find_row(self, arg):
        return to_str(super(cbf_wrapper, self).find_row(to_bytes(arg)))

    def set_datablockname(self, arg):
        return super(cbf_wrapper, self).set_datablockname(to_bytes(arg))

    def column_name(self):
        return to_str(super(cbf_wrapper, self).column_name())

    def get_value(self):
        return to_str(super(cbf_wrapper, self).get_value())

    def set_value(self, val):
        super(cbf_wrapper, self).set_value(to_bytes(val))

    def category_name(self):
        return to_str(super(cbf_wrapper, self).category_name())

    def get_typeofvalue(self):
        return to_str(super(cbf_wrapper, self).get_typeofvalue())

    def get_axis_type(self, axis_id):
        return to_str(super(cbf_wrapper, self).get_axis_type(to_bytes(axis_id)))

    def get_axis_offset(self, axis_id):
        return super(cbf_wrapper, self).get_axis_offset(to_bytes(axis_id))

    def get_axis_vector(self, axis_id):
        return super(cbf_wrapper, self).get_axis_vector(to_bytes(axis_id))

    def get_axis_setting(self, axis_id):
        return super(cbf_wrapper, self).get_axis_setting(to_bytes(axis_id))

    def get_axis_depends_on(self, axis_id):
        return to_str(super(cbf_wrapper, self).get_axis_depends_on(to_bytes(axis_id)))

    def get_axis_equipment_component(self, axis_id):
        return to_str(
            super(cbf_wrapper, self).get_axis_equipment_component(to_bytes(axis_id))
        )

    def set_axis_setting(self, axis_id, start, increment):
        super(cbf_wrapper, self).set_axis_setting(to_bytes(axis_id), start, increment)

    def construct_detector(self, element_number):
        return cbf_detector_wrapper(
            super(cbf_wrapper, self).construct_detector(element_number)
        )


class cbf_detector_wrapper(object):
    def __init__(self, cbf_detector):
        self._cbf_detector = cbf_detector

    def __getattr__(self, name):
        return getattr(self._cbf_detector, name)

    def __setattr__(self, name, value):
        if name == "_cbf_detector":
            super(cbf_detector_wrapper, self).__setattr__(name, value)
        else:
            setattr(self._cbf_detector, name, value)

    def get_detector_surface_axes(self, index):
        return to_str(self._cbf_detector.get_detector_surface_axes(index))

    def __swig_destroy__(self, cbf_detector):
        try:
            self._cbf_detector.__swig_destroy__(cbf_detector)
        except TypeError:
            self._cbf_detector.__swig_destroy__(cbf_detector._cbf_detector)
