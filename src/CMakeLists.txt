add_executable(umat_test main.cpp)

target_include_directories(umat_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(umat_test PRIVATE umat Model MCC SMCC)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/umat)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Model)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/models)