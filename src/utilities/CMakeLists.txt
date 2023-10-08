add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Integrator)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Intersection)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Math)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Tensor)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Types)

add_library(utilities INTERFACE)
target_include_directories(utilities INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(utilities INTERFACE plog Eigen3::Eigen Integrator Intersection Math Tensor Types)