# CMake requirements.
cmake_minimum_required(VERSION 3.5.0)
project(geomat
    VERSION 0.1.0
    LANGUAGES C CXX Fortran
)
include(CTest)
enable_testing()

# Settings.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wuninitialized -fsanitize=undefined,address -g")

# Verify CXX-Fortran interface compatibility.
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

# Add subdirectories of documentation and source code.
if(${BUILD_DOCS})
    add_subdirectory(${PROJECT_SOURCE_DIR}/docs)
endif()
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# Create Abaqus C++ interface integration test executable.
set(TARGET umat_cpp_test)
add_executable(${TARGET} tests/umat_integration.cpp)
target_link_libraries(${TARGET} PUBLIC umat)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create Abaqus Fortran interface integration test executable.
# Links against the standardU shared library to test the full Fortran→C++ path.
set(TARGET umat_fortran_test)
add_executable(${TARGET} tests/umat_integration.f90)
target_link_libraries(${TARGET} PUBLIC standardU)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)
set_target_properties(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_umat_fortran_test)

# Create static library for use with Abaqus, Plaxis etc...
add_library(geomat STATIC)
target_link_libraries(geomat PUBLIC umat)
target_compile_definitions(geomat PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Abaqus UMAT shared library (libstandardU.so / standardU.dll).
add_library(standardU SHARED
    src/umat/src/umat.f90
    src/umat/src/umat_cpp_interface.cpp
)
target_include_directories(standardU PRIVATE src/umat/src)
target_link_libraries(standardU PRIVATE abstract models utilities)
target_compile_definitions(standardU PRIVATE EIGEN_DONT_PARALLELIZE=1)
set_target_properties(standardU PROPERTIES
    OUTPUT_NAME "standardU"
    PREFIX "lib"
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_standardU
)

# Plaxis UDSM shared library (libusrmod.so / usrmod.dll).
add_library(usrmod SHARED
    src/udsm/src/udsm.f90
    src/udsm/src/udsm_cpp_interface.cpp
    src/udsm/src/udsm_metadata.cpp
)
target_include_directories(usrmod PRIVATE src/udsm/src)
target_link_libraries(usrmod PRIVATE abstract models utilities)
target_compile_definitions(usrmod PRIVATE EIGEN_DONT_PARALLELIZE=1)
set_target_properties(usrmod PROPERTIES
    OUTPUT_NAME "usrmod"
    PREFIX "lib"
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_usrmod
)

# UDSM C++ interface integration test executable.
set(TARGET udsm_cpp_test)
add_executable(${TARGET} tests/udsm_integration.cpp)
target_include_directories(${TARGET} PRIVATE src/udsm/src)
target_link_libraries(${TARGET} PUBLIC abstract models utilities usrmod)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create Python modules.
if(${pybind11_FOUND})
    pybind11_add_module(pyabstract MODULE src/geomat/pyabstract.cpp)
    target_link_libraries(pyabstract PUBLIC abstract utilities)
    target_compile_definitions(pyabstract PUBLIC EIGEN_DONT_PARALLELIZE=1)
    set_target_properties(pyabstract PROPERTIES OUTPUT_NAME "abstract")
    install(TARGETS pyabstract LIBRARY DESTINATION geomat)

    pybind11_add_module(pymodels MODULE src/geomat/pymodels.cpp)
    target_link_libraries(pymodels PUBLIC abstract models)
    target_compile_definitions(pymodels PUBLIC EIGEN_DONT_PARALLELIZE=1)
    set_target_properties(pymodels PROPERTIES OUTPUT_NAME "models")
    install(TARGETS pymodels LIBRARY DESTINATION geomat)

    pybind11_add_module(pyutilities MODULE src/geomat/pyutilities.cpp)
    target_link_libraries(pyutilities PUBLIC utilities)
    target_compile_definitions(pyutilities PUBLIC EIGEN_DONT_PARALLELIZE=1)
    set_target_properties(pyutilities PROPERTIES OUTPUT_NAME "utilities")
    install(TARGETS pyutilities LIBRARY DESTINATION geomat)
endif()

# Tests.
find_package(Catch2 QUIET)
if(TARGET Catch2::Catch2)
    message("Catch2 has been found. Tests will be generated.") 
    include(Catch)

    # Build tests.
    set(TARGET tests)
    add_executable(${TARGET} tests/tests.cpp)
    target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_link_libraries(tests Catch2::Catch2WithMain abstract models utilities umat usrmod)
    target_include_directories(tests PRIVATE src/udsm/src)
    catch_discover_tests(tests)
else()
    message("Catch2 not found hence no tests generated.")
endif()