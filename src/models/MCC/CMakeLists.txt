add_library(MCC INTERFACE)
target_sources(MCC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/MCC.cpp)
target_include_directories(MCC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(MCC INTERFACE Elastoplastic)