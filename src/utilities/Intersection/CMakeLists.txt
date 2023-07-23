add_library(Intersection INTERFACE)
target_sources(Intersection INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Intersection.cpp)
target_include_directories(Intersection INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Intersection INTERFACE plog Model Tensor Types)