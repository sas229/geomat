add_library(Integrators INTERFACE)
target_sources(Integrators INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Explicit/src/Explicit.cpp)
target_include_directories(Integrators INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/
    ${CMAKE_CURRENT_LIST_DIR}/src/Explicit/src
)
target_link_libraries(Integrators INTERFACE plog Types)