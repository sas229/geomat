cmake_minimum_required(VERSION 3.0.0)
project(umat_test VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DEIGEN_DONT_PARALLELIZE")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
message(${CMAKE_CXX_FLAGS})

add_subdirectory(${PROJECT_SOURCE_DIR}/docs)
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)