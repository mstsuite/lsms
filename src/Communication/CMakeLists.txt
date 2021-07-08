# Sources ----------------------------------------------------------------------
set(SRC
    distributeAtoms.cpp
    LSMSCommunication.cpp
    REWLCommunication.cpp
)

# Target -----------------------------------------------------------------------
add_library(Communication STATIC ${SRC})

# Libraries --------------------------------------------------------------------
target_link_libraries(Communication
    PUBLIC
        MPI::MPI_CXX
)

target_include_directories(Communication PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)