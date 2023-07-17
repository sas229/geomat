add_subdirectory(Explicit)

add_library(Integrator INTERFACE)
target_sources(Integrator INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/Integrator.cpp
)
target_include_directories(Integrator INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/
)
target_link_libraries(Integrator INTERFACE plog Types ModifiedEuler)