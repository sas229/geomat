# CMake requirements.
cmake_minimum_required(VERSION 3.0.0)
project(umat_test VERSION 0.1.0)
include(CTest)
enable_testing()

# Settings.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Add subdirectories of documentation and source code.
if (${BUILD_DOCS})
    add_subdirectory(${PROJECT_SOURCE_DIR}/docs)
endif()
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# Create test executable.
set(TARGET umat_test)
add_executable(${TARGET} src/main.cpp)
target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(${TARGET} PUBLIC umat Model C2MC LinearElastic MCC SMCC)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create library for usage with Abaqus.
set(TARGET umat_library)
add_library(${TARGET} SHARED src/umat/src/umat.cpp)
target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/umat)
target_link_libraries(${TARGET} PUBLIC Model C2MC LinearElastic MCC SMCC)
target_compile_definitions(${TARGET} PUBLIC EIGEN_DONT_PARALLELIZE=1)

# Create Python module.
set(TARGET models)
pybind11_add_module(${TARGET} src/geomat/models.cpp)
target_link_libraries(${TARGET} PUBLIC Types Model Elastic Elastoplastic C2MC LinearElastic MCC SMCC Eigen3::Eigen)
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