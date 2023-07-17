add_library(RKDP INTERFACE)
target_sources(RKDP INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/RKDP.cpp
)
target_include_directories(RKDP INTERFACE 
    ${CMAKE_CURRENT_LIST_DIR}/src/
)
target_link_libraries(RKDP INTERFACE plog Types Explicit)