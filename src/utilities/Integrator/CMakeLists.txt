add_library(Integrators INTERFACE)
target_sources(Integrators INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/Integrator.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/Explicit/ModifiedEuler/src/ModifiedEuler.cpp
)
target_include_directories(Integrators INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/
    ${CMAKE_CURRENT_LIST_DIR}/src/Explicit/ModifiedEuler/src/
)
target_link_libraries(Integrators INTERFACE plog Types)