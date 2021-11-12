# Distributes under BSD licence

#.rst:
# FindCCTBX
# ---------
#
# Find an existing CCTBX distribution and sets up targets for using it
#
# The CCTBX distribution is found by, in order:
#
# 1. Reading the ``CCTBX_BUILD_DIR`` cache variable, if set
# 2. Reading the ``LIBTBX_BUILD`` environment variable
# 3. Using python to ``import libtbx``
#
# Components
# ^^^^^^^^^^
#
# Any ``COMPONENTS`` passed to ``find_package`` will cause this module
# to look for a libtbx-distribution module of the same name. This will
# then be exposed as the target ``CCTBX::<module_name>`` which will set
# up the include, library paths associated with that module.
#
# Other libtbx module dependencies are not currently transitively
# followed, so will need to be specified manually.
#
# Any known shared libraries for a particular tbx-module will have
# targets added themselves. For example, cctbx builds a ``libcctbx``
# shared library. If the cctbx component is requested, this will be made
# available for linking with the target name ``CCTBX::cctbx::cctbx``.
#
# The database for recognising these shared library targets can be found
# in ``module_libraries.json`` in the folder above this FindCCTBX module.

# If python isn't already included, pull it in here
if (NOT TARGET Python::Interpreter)
    find_package(Python COMPONENTS Interpreter REQUIRED)
endif()

# Find the location of the libtbx build directory - where libtbx_env is
function(_cctbx_determine_libtbx_build_dir)
    # Try and read it from the environment
    message(DEBUG "Looking for libtbx build dir via environment LIBTBX_BUILD")
    if (DEFINED ENV{LIBTBX_BUILD})
        get_filename_component(_ABS_BUILD "$ENV{LIBTBX_BUILD}" ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}")
        set(LIBTBX_ENV "${_ABS_BUILD}/libtbx_env")
        message(DEBUG "Checking ${LIBTBX_ENV}")
        if (EXISTS "${LIBTBX_ENV}")
            set(CCTBX_BUILD_DIR "${_ABS_BUILD}" CACHE FILEPATH "Location of CCTBX build directory")
            message(DEBUG "Found libtbx via environment LIBTBX_BUILD: ${_ABS_BUILD}")
            return()
        endif()
    endif()

    message(DEBUG "Looking for libtbx build dir via importing libtbx in python")
    execute_process(COMMAND ${Python_EXECUTABLE} -c "import libtbx.load_env; print(abs(libtbx.env.build_path))"
                    RESULT_VARIABLE _LOAD_ENV_RESULT
                    OUTPUT_VARIABLE _LOAD_LIBTBX_BUILD_DIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    if (NOT ${_LOAD_ENV_RESULT})
        # We found it via python import
        message(DEBUG "Got libtbx build path: ${CCTBX_BUILD_DIR}")
        set(CCTBX_BUILD_DIR "${_LOAD_LIBTBX_BUILD_DIR}" CACHE FILEPATH "Location of CCTBX build directory")
        return()
    endif()

endfunction()

# Read details for a single module out of libtbx_env and other info
function(_cctbx_read_module MODULE)
    cmake_path(SET _read_env_script NORMALIZE "${CMAKE_CURRENT_LIST_DIR}/../read_env.py")
    execute_process(
        COMMAND ${Python_EXECUTABLE} "${_read_env_script}" "${CCTBX_BUILD_DIR}/libtbx_env"
        OUTPUT_VARIABLE _env_json
        RESULT_VARIABLE _result)
    if (_result)
        message(SEND_ERROR "Failed to read environment file: ${CCTBX_BUILD_DIR}/libtbx_env")
        return()
    endif()
    # We now have a json representation of libtbx_env - extract the entry for this modile
    string(JSON _module_json ERROR_VARIABLE _error GET "${_env_json}" module_dict ${MODULE})
    if (NOT _module_json)
        set("CCTBX_${MODULE}_FOUND" FALSE PARENT_SCOPE)
        return()
    else()
        set("CCTBX_${MODULE}_FOUND" TRUE PARENT_SCOPE)
    endif()
    # Now, ensure a target exists for this module
    set(_target "CCTBX::${MODULE}")
    if (NOT TARGET CCTBX::${MODULE})
        add_library(${_target} INTERFACE IMPORTED)
    endif()

    # Get the dist-wide include dir
    string(JSON _include_paths GET "${_env_json}" include_path)
    string(JSON _lib_path GET "${_env_json}" lib_path)

    # Work out what paths need to be included for this module
    string(JSON _n_dist_paths LENGTH "${_module_json}" dist_paths)
    math(EXPR _n_dist_paths "${_n_dist_paths} - 1")
    # We need to account for:
    # Algorithm: if folder has an include/ subdir:
    #               use that
    #            else:
    #               use the path above
    # - this accounts for most cases, and it seems unlikely that we will
    #   be importing uniquely from modules that this doesn't cover. There
    #   was an exhaustive mapping list in tbx2cmake but the new installed
    #   layout confuses that somewhat.
    foreach(_n RANGE "${_n_dist_paths}")
        string(JSON _dist_path GET "${_module_json}" dist_paths ${_n})
        if(NOT _dist_path)
            continue()
        endif()
        list(APPEND _dist_paths "${_dist_path}/include")
        if (EXISTS "${_dist_path}/include")
            list(APPEND _include_paths "${_dist_path}/include")
        else()
            cmake_path(GET _dist_path PARENT_PATH _result)
            list(APPEND _include_paths "${_result}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _include_paths)
    list(REMOVE_DUPLICATES _dist_paths)

    target_include_directories(${_target} INTERFACE ${_include_paths})

    # Find if this module has a "base" library, and any sub-libraries
    file(READ "${CMAKE_CURRENT_LIST_DIR}/../module_libraries.json" _modules_libs)
    string(JSON _module_libs GET "${_modules_libs}" "${MODULE}")
    if (_module_libs)
        # We have actual libraries to import as imported libraries
        # iterate over every key: value in the object
        string(JSON _n_libs LENGTH "${_module_libs}")
        math(EXPR _n_libs "${_n_libs} - 1")
        foreach(_n RANGE ${_n_libs})
            string(JSON _name MEMBER "${_module_libs}" ${_n})
            string(JSON _libname GET "${_module_libs}" "${_name}")
            # Find this library
            find_library(_liblocation ${_libname} HINTS ${_lib_path})
            if(_name STREQUAL base)
                message(DEBUG "Module ${_target} has root library ${_libname}")
                if (NOT _liblocation)
                    message(WARNING "Libtbx module ${MODULE} has base library named lib${_libname} but cannot find it - the module may be misconfigured")
                else()
                    target_link_libraries(${_target} INTERFACE ${_liblocation})
                endif()
            else()
                message(DEBUG "Extra library ${_target}::${_name} = ${_libname}")
                add_library(${_target}::${_name} INTERFACE IMPORTED)
                target_link_libraries(${_target}::${_name} INTERFACE ${_target} ${_libname})
            endif()
        endforeach()
    endif()
endfunction()


if (NOT CCTBX_BUILD_DIR)
    _cctbx_determine_libtbx_build_dir()
endif()

if (CCTBX_BUILD_DIR)
    message(DEBUG "Using build dir ${CCTBX_BUILD_DIR}")

    foreach(_comp IN LISTS CCTBX_FIND_COMPONENTS)
        # Only try and find this component if we didn't already find it
        if(NOT TARGET CCTBX::${_comp})
            _cctbx_read_module(${_comp})
        else()
            # We already had a target for this module, so mark found
            set("CCTBX_${_comp}_FOUND" TRUE)
        endif()
    endforeach()
    unset(_comp)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CCTBX
    REQUIRED_VARS CCTBX_BUILD_DIR
    HANDLE_COMPONENTS
)
