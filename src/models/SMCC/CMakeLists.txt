add_library(SMCC INTERFACE)
target_sources(SMCC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/SMCC.cpp)
target_include_directories(SMCC INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(SMCC INTERFACE Elastoplastic)