# Find Sphinx uing custom cmake file.
find_package(Sphinx REQUIRED)

# Set source and output directories.
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/xml)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx_html)

# Build documentation.
add_custom_target(Sphinx ALL
    COMMAND ${SPHINX_EXECUTABLE} -b html
    # Tell Breathe where to find the Doxygen output.
    -Dbreathe_projects.geomat=${DOXYGEN_OUTPUT_DIR}
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx...")