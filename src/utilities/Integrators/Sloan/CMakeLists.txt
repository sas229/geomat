add_library(Sloan INTERFACE)
target_sources(Sloan INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Sloan.cpp)
target_include_directories(Sloan INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Sloan INTERFACE plog Types)