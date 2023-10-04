
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/C2MC)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/EMC)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/LinearElastic)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/MCC)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/SMCC)

add_library(models INTERFACE)
target_include_directories(models INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(models INTERFACE C2MC EMC LinearElastic MCC SMCC)