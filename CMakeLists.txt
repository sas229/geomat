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

# Verify CXX Fortran interface compatibility.
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

# Add subdirectories of documentation and source code.
if (${BUILD_DOCS})
    add_subdirectory(${PROJECT_SOURCE_DIR}/docs)
endif()
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# Create Fortran interface test executable.
set(TARGET umat_fortran_test)
add_executable(${TARGET} tests/umat_integration.f90)
target_link_libraries(${TARGET} PUBLIC umat Model Models)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create C++ test executable.
set(TARGET umat_cpp_test)
add_executable(${TARGET} tests/umat_integration.cpp)
target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(${TARGET} PUBLIC umat Model Models)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create Python module.
set(TARGET models)
pybind11_add_module(${TARGET} src/geomat/models.cpp)
target_link_libraries(${TARGET} PUBLIC Types Model Elastic Elastoplastic Models Eigen3::Eigen)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)
install(TARGETS models DESTINATION ./)

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