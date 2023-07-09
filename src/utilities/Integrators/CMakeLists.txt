add_subdirectory(Sloan)

add_library(Integrators INTERFACE)
target_include_directories(Integrators INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Integrators INTERFACE Sloan)