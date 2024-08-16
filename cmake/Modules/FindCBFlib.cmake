
# Distributes under BSD licence

#.rst:
# FindCBFlib
# ----------
#
# Find CBFlib and create an imported target for it.
#
# If they cannot be found in the standard search, you may need to add
# hints of:
#
# ``CBFLIB_CBF_LIB``
#   The location of the CBFlib library file
# ``CBFLIB_CBF_H``
#   The path to the cbflib.h include file
#
# The imported target ``CBFlib::cbf`` will be created.
#
# Img
# ^^^
#
# If it can be found, optionally via the ``CBFLIB_IMG_H`` and
# ``CBFLIB_CBF_H`` cache variables, the target ``CBFlib::img`` will be
# added for linking to libimg. This can be specified as a required
# component.

find_file(CBFLIB_CBF_H cbf.h PATH_SUFFIXES cbflib)
find_library(CBFLIB_CBF_LIB cbf)

find_file(CBFLIB_IMG_H img.h PATH_SUFFIXES cbflib)
find_library(CBFLIB_IMG_LIB cbf)

if (CBFLIB_CBF_H AND CBFLIB_CBF_LIB)
    set(CBFLIB_CBF_FOUND TRUE)
    set(CBFLIB_FOUND TRUE)
endif()

if (CBFLIB_IMG_H AND CBFLIB_IMG_LIB)
    set(CBFLIB_IMG_FOUND TRUE)
endif()

if (CBFLIB_CBF_FOUND AND NOT TARGET CBFlib::cbf)
    add_library(CBFlib::cbf UNKNOWN IMPORTED)
    cmake_path(GET CBFLIB_CBF_H PARENT_PATH CBFLIB_CBF_DIR)
    set_target_properties(CBFlib::cbf PROPERTIES IMPORTED_LOCATION ${CBFLIB_CBF_LIB})
    target_include_directories(CBFlib::cbf INTERFACE "${CBFLIB_CBF_DIR}")
endif()

if (CBFLIB_IMG_FOUND AND NOT TARGET CBFlib::img)
    add_library(CBFlib::img UNKNOWN IMPORTED)
    set_target_properties(CBFlib::img PROPERTIES IMPORTED_LOCATION ${CBFLIB_IMG_LIB})
    cmake_path(GET CBFLIB_IMG_H PARENT_PATH CBFLIB_IMG_DIR)
    target_include_directories(CBFlib::img INTERFACE "${CBFLIB_IMG_DIR}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBFlib
    REQUIRED_VARS CBFLIB_CBF_H CBFLIB_CBF_LIB
    HANDLE_COMPONENTS
)
