add_library(Tensor INTERFACE)
target_sources(Tensor INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Tensor.cpp)
target_include_directories(Tensor INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Tensor INTERFACE plog Eigen3::Eigen Types)