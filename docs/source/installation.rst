Installation
============

Prerequisites
-------------

- CMake 3.5 or later
- C/C++ compiler (GCC, Clang, or MSVC)
- Fortran compiler (gfortran or Intel Fortran) — required for the Abaqus and Plaxis shared libraries
- Python 3.6+ with development headers (optional, for Python bindings)

Building
--------

.. code-block:: console

    git clone https://github.com/sas229/geomat.git
    cd geomat
    git submodule init src/external
    git submodule update --init --recursive
    mkdir build && cd build
    cmake ..
    make

Build Targets
-------------

The following targets are available:

- ``standardU`` — Abaqus UMAT shared library (``libstandardU.so`` / ``standardU.dll``)
- ``usrmod`` — Plaxis UDSM shared library (``libusrmod.so`` / ``usrmod.dll``)
- ``geomat`` — Static library
- ``pyabstract``, ``pymodels``, ``pyutilities`` — Python modules (requires pybind11)
- ``tests`` — Catch2 unit tests
- ``umat_cpp_test`` — UMAT C++ integration test
- ``umat_fortran_test`` — UMAT Fortran integration test
- ``udsm_cpp_test`` — UDSM C++ integration test

To build a specific target:

.. code-block:: console

    make standardU

CMake Options
-------------

- ``BUILD_DOCS`` — Build Sphinx/Doxygen documentation (requires Sphinx and Breathe)
