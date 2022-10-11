
#
#
#
add_executable(LIZandVoronoiTest)
set_target_properties(LIZandVoronoiTest PROPERTIES OUTPUT_NAME "LIZandVoronoiTest")

set_target_properties(LIZandVoronoiTest
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

target_link_libraries(LIZandVoronoiTest
        PUBLIC lsmscore)

#
#
target_sources(LIZandVoronoiTest
        PUBLIC
        LIZandVoronoiTest.cpp)

