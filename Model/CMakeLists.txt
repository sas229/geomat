add_library(Model INTERFACE)
find_package(Eigen3 3.3 REQUIRED)
target_sources(Model INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Model.cpp)
target_include_directories(Model INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)