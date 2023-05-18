import sys
import os
import libtbx.load_env
from libtbx.env_config import get_boost_library_with_python_version
from pathlib import Path

try:
    # Check if we have a "modernized" pycbf
    import pycbf

    if pycbf.__version__.startswith("CBFlib"):
        # Old pycbf
        raise ImportError
    build_cbf_bindings = False
except ImportError:
    # We have to build CBFlib bindings ourselves
    build_cbf_bindings = True

Import("env_etc")

env_etc.dxtbx_dist = libtbx.env.dist_path("dxtbx")
env_etc.dxtbx_include = os.path.join(env_etc.dxtbx_dist, "src")
env_etc.dxtbx_includes = []
env_etc.dxtbx_common_includes = [
    env_etc.base_include,
    env_etc.libtbx_include,
    env_etc.scitbx_include,
    env_etc.boost_adaptbx_include,
    env_etc.boost_include,
    env_etc.dxtbx_include,
]

Import("env_no_includes_boost_python_ext")
env = env_no_includes_boost_python_ext.Clone()

if libtbx.env.build_options.use_conda:
    boost_python = get_boost_library_with_python_version(
        "boost_python", env_etc.conda_libpath
    )
    env.Append(LIBPATH=env_etc.conda_libpath)
    env_etc.dxtbx_common_includes.extend(env_etc.conda_cpppath)
else:
    boost_python = "boost_python"

env_etc.dxtbx_libs = ["tiff", boost_python]
env_etc.dxtbx_lib_paths = [
    env_etc.base_lib,
    env_etc.libtbx_lib,
    os.path.join(sys.prefix, "lib"),
]

env_etc.dxtbx_hdf5_libs = ["hdf5"]
env_etc.dxtbx_hdf5_lib_paths = []

if sys.platform == "win32" and env_etc.compiler == "win32_cl":
    if libtbx.env.build_options.use_conda:
        env_etc.dxtbx_hdf5_libs = ["hdf5"]
        env_etc.cppdefines = {"H5_BUILT_AS_DYNAMIC_LIB": 1}
        env_etc.dxtbx_hdf5_libs.append("zlib")
        env_etc.dxtbx_includes.extend(env_etc.conda_cpppath)
        env_etc.dxtbx_lib_paths.extend(env_etc.conda_libpath)
    else:
        env_etc.dxtbx_lib_paths = [
            env_etc.libpath_python,
            env_etc.libtbx_lib,
            libtbx.env.under_base("libtiff"),
        ]
        env_etc.dxtbx_hdf5_lib_paths = [
            libtbx.env.under_base(os.path.join("HDF5-1.8.16", "lib"))
        ]
        env_etc.dxtbx_includes.append(
            libtbx.env.under_base(os.path.join("HDF5-1.8.16", "include"))
        )
        env_etc.dxtbx_includes.append(
            os.path.join(env_etc.cctbx_include, "msvc9.0_include")
        )
        env_etc.dxtbx_includes.append(libtbx.env.under_base("libtiff"))
        env_etc.dxtbx_libs = ["libtiff", "boost_python"]

if build_cbf_bindings:
    env_etc.dxtbx_common_includes.extend(env_etc.cbflib_common_includes)
    env_etc.dxtbx_libs.append("cbf")

# for the hdf5.h file - look at where Python is coming from unless is OS X
# framework build... messy but appears to work on Linux and OS X
include_root = os.path.split(env_etc.python_include)[0]
if "Python.framework" in include_root:
    include_root = os.path.join(include_root.split("Python.framework")[0], "include")
if os.path.exists(os.path.join(include_root, "hdf5.h")):
    env_etc.dxtbx_includes.append(include_root)
else:
    # check for PSDM installation. Example:
    # /reg/g/psdm/sw/external/hdf5/1.8.6/x86_64-rhel5-gcc41-opt/include
    psdm_hdf5_path = os.path.join(
        os.environ.get("SIT_ROOT", ""),
        "sw",
        "external",
        "hdf5",
        "1.8.6",
        os.environ.get("SIT_ARCH", ""),
        "include",
    )
    if os.path.exists(psdm_hdf5_path):
        env_etc.dxtbx_common_includes.append(psdm_hdf5_path)
    psdm_hdf5_path = os.path.join(
        os.environ.get("SIT_ROOT", ""),
        "sw",
        "external",
        "hdf5",
        "1.8.6",
        os.environ.get("SIT_ARCH", ""),
        "lib",
    )
    if os.path.exists(psdm_hdf5_path):
        env_etc.dxtbx_hdf5_lib_paths.append(psdm_hdf5_path)

if not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include"):
    Import("env_no_includes_boost_python_ext")
    env = env_no_includes_boost_python_ext.Clone()

    # Don't surface warnings from system or cctbx_project headers
    system_includes = [x for x in env_etc.conda_cpppath if x]
    system_includes.append(str(Path(env_etc.scitbx_dist).parent))
    env.Append(CXXFLAGS=[f"-isystem{x}" for x in system_includes])
    env.Append(SHCXXFLAGS=[f"-isystem{x}" for x in system_includes])

    env_etc.enable_more_warnings(env=env)
    env_etc.include_registry.append(
        env=env,
        paths=env_etc.dxtbx_includes
        + env_etc.dxtbx_common_includes
        + [env_etc.python_include],
    )

    env.Append(
        LIBS=env_etc.libm + ["cctbx", "scitbx_boost_python"] + env_etc.dxtbx_libs,
        LIBPATH=env_etc.dxtbx_lib_paths + env_etc.dxtbx_hdf5_lib_paths,
    )

    # Fix the build environment so that it doesn't break on modern C++
    for path in list(env["CPPPATH"]):
        if "msvc9.0_include" in path:
            env["CPPPATH"].remove(path)

    if env_etc.clang_version:
        wd = ["-Wno-unused-function"]
        env.Append(CCFLAGS=wd)

    if build_cbf_bindings:
        env.Append(CPPDEFINES="BUILD_CBF")

    if hasattr(env_etc, "cppdefines"):
        env.Append(CPPDEFINES=env_etc.cppdefines)

    env.SharedLibrary(
        target="#lib/dxtbx_ext",
        source=[
            "src/dxtbx/boost_python/ext.cpp",
            "src/dxtbx/boost_python/compression.cc",
        ],
        LIBS=env_etc.libs_python + env_etc.libm + env_etc.dxtbx_libs,
    )

    nexus = env.SharedLibrary(
        target="#/lib/dxtbx_format_nexus_ext",
        source=["src/dxtbx/format/boost_python/nexus_ext.cc"],
        LIBS=env_etc.libs_python
        + env_etc.libm
        + env_etc.dxtbx_libs
        + env_etc.dxtbx_hdf5_libs,
    )

    imageset = env.SharedLibrary(
        target="#/lib/dxtbx_imageset_ext",
        source=["src/dxtbx/boost_python/imageset_ext.cc"],
        LIBS=env_etc.libs_python
        + env_etc.libm
        + env_etc.dxtbx_libs
        + env_etc.dxtbx_hdf5_libs,
    )

    dxtbx_format_image_ext_sources = [
        "src/dxtbx/format/boost_python/image_ext.cc",
    ]
    if build_cbf_bindings:
        dxtbx_format_image_ext_sources.append(
            "src/dxtbx/format/boost_python/cbf_read_buffer.cpp"
        )
    image = env.SharedLibrary(
        target="#/lib/dxtbx_format_image_ext",
        source=dxtbx_format_image_ext_sources,
        LIBS=env_etc.libs_python
        + env_etc.libm
        + env_etc.dxtbx_libs
        + env_etc.dxtbx_hdf5_libs,
    )

    model = env.SharedLibrary(
        target="#/lib/dxtbx_model_ext",
        source=[
            "src/dxtbx/model/boost_python/beam.cc",
            "src/dxtbx/model/boost_python/spectrum.cc",
            "src/dxtbx/model/boost_python/goniometer.cc",
            "src/dxtbx/model/boost_python/kappa_goniometer.cc",
            "src/dxtbx/model/boost_python/multi_axis_goniometer.cc",
            "src/dxtbx/model/boost_python/panel.cc",
            "src/dxtbx/model/boost_python/detector.cc",
            "src/dxtbx/model/boost_python/scan.cc",
            "src/dxtbx/model/boost_python/scan_helpers.cc",
            "src/dxtbx/model/boost_python/crystal.cc",
            "src/dxtbx/model/boost_python/parallax_correction.cc",
            "src/dxtbx/model/boost_python/pixel_to_millimeter.cc",
            "src/dxtbx/model/boost_python/experiment.cc",
            "src/dxtbx/model/boost_python/experiment_list.cc",
            "src/dxtbx/model/boost_python/model_ext.cc",
        ],
        LIBS=env_etc.libs_python + env_etc.libm + env_etc.dxtbx_libs + env["LIBS"],
    )

    flumpy_flags = []
    # cl.exe treats symbols this way by default
    if env_etc.compiler != "win32_cl":
        flumpy_flags = ["-fvisibility=hidden"]

    env.SharedLibrary(
        target="#/lib/dxtbx_flumpy",
        source=[
            "src/dxtbx/boost_python/flumpy.cc",
        ],
        LIBS=env_etc.libs_python + env_etc.libm + env_etc.dxtbx_libs + env["LIBS"],
        CPPFLAGS=flumpy_flags,
    )

    env.SConscript("src/dxtbx/masking/SConscript", exports={"env": env})
