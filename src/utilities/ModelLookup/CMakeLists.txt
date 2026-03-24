add_library(ModelLookup INTERFACE)
target_sources(ModelLookup INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/ModelLookup.cpp)
target_include_directories(ModelLookup INTERFACE ${CMAKE_CURRENT_LIST_DIR}/src/)
target_link_libraries(ModelLookup INTERFACE abstract models)
