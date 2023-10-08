# CMake requirements.
cmake_minimum_required(VERSION 3.5.0)
project(geomat 
    VERSION 0.1.0
    LANGUAGES C CXX
)
include(CTest)
enable_testing()

# Settings.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Add subdirectories of documentation and source code.
if(${BUILD_DOCS})
    add_subdirectory(${PROJECT_SOURCE_DIR}/docs)
endif()
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# if(NOT WIN32)
#     # Verify CXX Fortran interface compatibility.
#     enable_language(Fortran)
#     include(FortranCInterface)
#     FortranCInterface_VERIFY(CXX)

#     # Create Abaqus Fortran interface integration test executable.
#     set(TARGET umat_fortran_test)
#     add_executable(${TARGET} tests/umat_integration.f90)
#     target_link_libraries(${TARGET} PUBLIC umat)
#     target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)
# endif()

# Create Abaqus C++ interface integration test executable.
set(TARGET umat_cpp_test)
add_executable(${TARGET} tests/umat_integration.cpp)
target_link_libraries(${TARGET} PUBLIC umat)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create static library for use with Abaqus, Plaxis etc...
add_library(geomat STATIC)
target_link_libraries(geomat PUBLIC umat)
target_compile_definitions(geomat PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create Python modules.
if(${pybind11_FOUND})
    pybind11_add_module(pyabstract MODULE src/geomat/pyabstract.cpp)
    target_link_libraries(pyabstract PUBLIC abstract utilities)
    target_compile_definitions(pyabstract PUBLIC EIGEN_DONT_PARALLELIZE=1)
    set_target_properties(pyabstract PROPERTIES OUTPUT_NAME "abstract")
    install(TARGETS pyabstract LIBRARY DESTINATION ./)

    pybind11_add_module(pymodels MODULE src/geomat/pymodels.cpp)
    target_link_libraries(pymodels PUBLIC abstract models)
    target_compile_definitions(pymodels PUBLIC EIGEN_DONT_PARALLELIZE=1)
    set_target_properties(pymodels PROPERTIES OUTPUT_NAME "models")
    install(TARGETS pymodels LIBRARY DESTINATION ./)

    pybind11_add_module(pyutilities MODULE src/geomat/pyutilities.cpp)
    target_link_libraries(pyutilities PUBLIC utilities)
    target_compile_definitions(pyutilities PUBLIC EIGEN_DONT_PARALLELIZE=1)
    set_target_properties(pyutilities PROPERTIES OUTPUT_NAME "utilities")
    install(TARGETS pyutilities LIBRARY DESTINATION ./)
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
    target_link_libraries(tests Catch2::Catch2WithMain MCC)
    catch_discover_tests(tests)
else()
    message("Catch2 not found hence no tests generated.")
endif()