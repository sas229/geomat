
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/C2MC)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/EMC)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/LinearElastic)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/MCC)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/SMCC)

add_library(Models INTERFACE)
target_link_libraries(Models INTERFACE C2MC EMC LinearElastic MCC SMCC)