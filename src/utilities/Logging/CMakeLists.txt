add_library(Logging INTERFACE)
target_sources(Logging INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Logging.cpp)
target_include_directories(Logging INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Logging INTERFACE plog)