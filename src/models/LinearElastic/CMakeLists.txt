add_library(LinearElastic INTERFACE)
target_sources(LinearElastic INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/LinearElastic.cpp)
target_include_directories(LinearElastic INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(LinearElastic INTERFACE Elastic)