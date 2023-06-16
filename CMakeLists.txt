# CMake requirements.
cmake_minimum_required(VERSION 3.0.0)
project(umat_test VERSION 0.1.0)
include(CTest)
enable_testing()

# Settings.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DEIGEN_DONT_PARALLELIZE")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
message("Compiler flags:" ${CMAKE_CXX_FLAGS})

# Add subdirectories of documentation and source code.
# add_subdirectory(${PROJECT_SOURCE_DIR}/docs)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

# # Create test executable.
# add_executable(umat_test src/main.cpp)
# target_include_directories(umat_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
# target_link_libraries(umat_test PUBLIC umat Model LinearElastic MCC SMCC)

# # Create library for usage with Abaqus.
# add_library(umat_library SHARED src/umat/src/umat.cpp)
# target_include_directories(umat_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/umat)
# target_link_libraries(umat_library PUBLIC Model LinearElastic MCC SMCC)

# Create Python module.
pybind11_add_module(models src/geomat/models.cpp)
target_link_libraries(models PUBLIC Types Model Elastic Elastoplastic LinearElastic MCC SMCC Eigen3::Eigen)
install(TARGETS models DESTINATION ./)

# Tests.
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)