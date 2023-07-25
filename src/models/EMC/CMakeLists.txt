add_library(EMC INTERFACE)
target_sources(EMC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/EMC.cpp)
target_include_directories(EMC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(EMC INTERFACE Elastoplastic)