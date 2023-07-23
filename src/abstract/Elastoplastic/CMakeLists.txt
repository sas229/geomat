add_library(Elastoplastic INTERFACE)
target_sources(Elastoplastic INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Elastoplastic.cpp)
target_include_directories(Elastoplastic INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Elastoplastic INTERFACE Elastic)