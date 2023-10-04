add_library(umat INTERFACE)
target_sources(umat INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/umat_cpp_interface.cpp)
target_include_directories(umat INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(umat INTERFACE abstract models utilities)