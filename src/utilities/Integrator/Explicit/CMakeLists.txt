add_library(Explicit INTERFACE)
target_sources(Explicit INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/Explicit.cpp
)
target_include_directories(Explicit INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/
)
target_link_libraries(Explicit INTERFACE plog Types)

add_subdirectory(ModifiedEuler)
add_subdirectory(RKDP)