set(SRC
    getcor_c.f
    calculateCoreStates.cpp
    corslv_c.f
    deepst_c.f
    semcst_c.f
    inws_c.f
    outws_c.f
    inwhnk_c.f
    invals_c.f
    richnk_c.f
)

# Target -----------------------------------------------------------------------
add_library(Core STATIC ${SRC})

# Libraries --------------------------------------------------------------------
#target_link_libraries(Core
#    PUBLIC
#        ${BLAS_LIBRARIES}
#        ${LAPACK_LIBRARIES}
#)

target_include_directories(Core PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)