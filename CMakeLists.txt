cmake_minimum_required(VERSION 3.0.0)
project(umat_test VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DEIGEN_DONT_PARALLELIZE")
message(${CMAKE_CXX_FLAGS})

add_subdirectory(${PROJECT_SOURCE_DIR}/umat)
add_subdirectory(${PROJECT_SOURCE_DIR}/Model)
add_subdirectory(${PROJECT_SOURCE_DIR}/models)

add_executable(umat_test main.cpp)

target_include_directories(umat_test PRIVATE ${PROJECT_SOURCE_DIR})

target_link_libraries(umat_test PRIVATE umat Model MCC SMCC)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
