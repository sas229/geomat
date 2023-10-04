add_library(Model INTERFACE)
target_sources(Model INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Model.cpp)
target_include_directories(Model INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Model INTERFACE plog Eigen3::Eigen utilities)