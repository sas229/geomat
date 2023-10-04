add_library(Elastic INTERFACE)
target_sources(Elastic INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/Elastic.cpp)
target_include_directories(Elastic INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(Elastic INTERFACE Model)