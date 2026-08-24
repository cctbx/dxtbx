#!/usr/bin/env bash

if [[ $CONDA_TOOLCHAIN_BUILD != $CONDA_TOOLCHAIN_HOST ]]; then
    # Conda does some swizzling when cross compiling, including moving
    # the site-packages folder to the build prefix. So let's just
    # manually add this to the compiler search path.
    CFLAGS="-isystem $BUILD_PREFIX/lib/python$PY_VER/site-packages $CXXFLAGS"
    CXXFLAGS="-isystem $BUILD_PREFIX/lib/python$PY_VER/site-packages $CXXFLAGS"
fi

mkdir _build
cd _build
cmake ${CMAKE_ARGS} ../dxtbx "-DCMAKE_INSTALL_PREFIX=$PREFIX" "-DPython_EXECUTABLE=$PYTHON" -GNinja
cmake --build .
cmake --install .
$PYTHON -mpip install ../dxtbx