
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

Python_add_library( dxtbx_ext
    MODULE
    boost_python/ext.cpp
    boost_python/compression.cc )
target_link_libraries(dxtbx_ext PRIVATE Boost::python CCTBX::scitbx)

Python_add_library( dxtbx_format_nexus_ext MODULE format/boost_python/nexus_ext.cc )
target_link_libraries( dxtbx_format_nexus_ext PRIVATE  Boost::python CCTBX::scitbx hdf5::hdf5 )

Python_add_library( dxtbx_imageset_ext MODULE boost_python/imageset_ext.cc )
target_link_libraries( dxtbx_imageset_ext PUBLIC Boost::python CCTBX::scitbx )

Python_add_library( dxtbx_format_image_ext
    MODULE
    format/boost_python/image_ext.cc
)
target_link_libraries( dxtbx_format_image_ext PUBLIC Boost::python CCTBX::scitbx )
# If we have CBFlib available, use it
if (TARGET CBFlib::cbf)
    target_link_libraries( dxtbx_format_image_ext PRIVATE CBFlib::cbf )
    target_compile_definitions( dxtbx_format_image_ext PRIVATE BUILD_CBF )
    target_sources( dxtbx_format_image_ext PRIVATE format/boost_python/cbf_read_buffer.cpp )
endif()

Python_add_library( dxtbx_model_ext
    MODULE
    model/boost_python/beam.cc
    model/boost_python/spectrum.cc
    model/boost_python/goniometer.cc
    model/boost_python/kappa_goniometer.cc
    model/boost_python/multi_axis_goniometer.cc
    model/boost_python/panel.cc
    model/boost_python/detector.cc
    model/boost_python/scan.cc
    model/boost_python/scan_helpers.cc
    model/boost_python/crystal.cc
    model/boost_python/parallax_correction.cc
    model/boost_python/pixel_to_millimeter.cc
    model/boost_python/experiment.cc
    model/boost_python/experiment_list.cc
    model/boost_python/model_ext.cc )
target_link_libraries(
    dxtbx_model_ext
    PUBLIC
    Boost::python
    CCTBX::scitbx
)
if(WIN32)
    # Technically required on other platforms, but already loaded by time
    # of use in dxtbx, and a bug with -install_name settings for these
    # prevents direct linking from working without fiddling around with
    # install_name_tool.
    target_link_libraries(
        dxtbx_model_ext
        PUBLIC
        CCTBX::scitbx::boost_python
        CCTBX::cctbx
    )
endif()

Python_add_library( dxtbx_masking_ext
    MODULE
    masking/boost_python/ext.cc )
target_link_libraries( dxtbx_masking_ext PUBLIC Boost::python CCTBX::scitbx )

pybind11_add_module(dxtbx_flumpy boost_python/flumpy.cc)
target_link_libraries(dxtbx_flumpy PUBLIC Boost::python CCTBX::scitbx )

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fvisibility=hidden HAS_VISIBILITY)
if(HAS_VISIBILITY)
    target_compile_options(dxtbx_flumpy PRIVATE -fvisibility=hidden)
endif()

install(
    TARGETS
        dxtbx_ext
        dxtbx_format_nexus_ext
        dxtbx_imageset_ext
        dxtbx_format_image_ext
        dxtbx_model_ext
        dxtbx_masking_ext
        dxtbx_flumpy
    DESTINATION "${Python_SITEARCH}"
)

