add_library(Checks INTERFACE)
target_sources(Checks INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Checks.cpp)
target_include_directories(Checks INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Checks INTERFACE plog)