cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(dxtbx)

# Add the included modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# General cmake environment configuration
include(SetDefaultBuildRelWithDebInfo) # Default builds to release with debug info
include(AlwaysColourCompilation) # Always show coloured compiler output
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Generate compile_commands.json
set(CMAKE_CXX_STANDARD 11)

find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(CCTBX COMPONENTS scitbx REQUIRED)
find_package(HDF5 REQUIRED)
find_package(pybind11 REQUIRED)
find_package(CBFlib)

# Find the boost::python library for this version of python
find_package(Boost COMPONENTS "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}" REQUIRED)

# Create Boost::python alias so we don't need to carry the python version around
if(NOT TARGET Boost::python )
    message("Adding Boost::python target")
    add_library(Boost::python INTERFACE IMPORTED)
    set_target_properties(Boost::python PROPERTIES
        INTERFACE_LINK_LIBRARIES  "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}" )
endif()

# Put the libraries into lib/ so that we can run this in-place in a TBX install
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/lib")

add_subdirectory(src/dxtbx)

