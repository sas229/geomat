add_library(Math INTERFACE)
target_sources(Math INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Math.cpp)
target_include_directories(Math INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Math INTERFACE plog Types)