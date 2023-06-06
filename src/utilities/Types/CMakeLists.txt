add_library(Types INTERFACE)
target_sources(Types INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Types.cpp)
target_include_directories(Types INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Types INTERFACE plog Eigen3::Eigen)