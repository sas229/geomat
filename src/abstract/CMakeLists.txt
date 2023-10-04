add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Model)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Elastic)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Elastoplastic)

add_library(abstract INTERFACE)
target_include_directories(abstract INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(abstract INTERFACE Model Elastic Elastoplastic)
