add_library(C2MC INTERFACE)
target_sources(C2MC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/C2MC.cpp)
target_include_directories(C2MC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(C2MC INTERFACE Checks Logging Types Elastoplastic)