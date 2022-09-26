cmake_minimum_required(VERSION 3.0.0)
project(umat_test VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DEIGEN_DONT_PARALLELIZE")
message(${CMAKE_CXX_FLAGS})

find_package (Eigen3 3.3 REQUIRED)
# find_package (Boost REQUIRED)

add_executable(umat_test)

target_sources(umat_test PRIVATE  main.cpp)

target_link_libraries(umat_test PUBLIC Eigen3::Eigen)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
