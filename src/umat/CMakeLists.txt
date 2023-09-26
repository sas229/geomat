add_library(umat INTERFACE)
target_sources(umat 
    INTERFACE 
        # ${CMAKE_CURRENT_LIST_DIR}/src/umat.f90
        # ${CMAKE_CURRENT_LIST_DIR}/src/abaqus_fortran_cpp_interface.f90
        ${CMAKE_CURRENT_LIST_DIR}/src/umat_cpp_interface.cpp
        # ${CMAKE_CURRENT_LIST_DIR}/src/umat.cpp
)
target_include_directories(umat INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(umat INTERFACE plog Eigen3::Eigen Types Model Models)

# add_library(umat INTERFACE)
# target_sources(umat INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/umat.cpp)
# target_include_directories(umat INTERFACE 
#     ${CMAKE_CURRENT_LIST_DIR}/src/)
# target_link_libraries(umat INTERFACE plog Eigen3::Eigen Types Model Models)