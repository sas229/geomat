add_library(ModifiedEuler INTERFACE)
target_sources(ModifiedEuler INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/ModifiedEuler.cpp
)
target_include_directories(ModifiedEuler INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/
)
target_link_libraries(ModifiedEuler INTERFACE plog Tensor Types Explicit)