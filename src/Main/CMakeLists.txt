# Source Lists -----------------------------------------------------------------
set(LIB_LSMS_SOURCES
    SystemParameters.cpp
    read_input.cpp
    PotentialIO.cpp
    buildLIZandCommLists.cpp
    energyContourIntegration.cpp
    solveSingleScatterers.cpp
    calculateDensities.cpp
    calculateChemPot.cpp
    checkConsistency.cpp
    lsmsClass.cpp
    calculateEvec.cpp
    initializeAtom.cpp
    mixing.cpp
    ReplicaExchangeWL.cpp
    AlloyBankIO.cpp
    rotateToGlobal.cpp
    write_restart.cpp
)

# Libraries --------------------------------------------------------------------
add_library(liblsms STATIC ${LIB_LSMS_SOURCES})

set(LIBS
        Communication
        RadialGrid
        MultipleScattering
        VORPOL
        Madelung
        Potential
        SingleSite
        Core
        LuaInterface
        Misc
        Accelerator
        TotalEnergy
        ${HDF5_LIBRARIES}
)

target_include_directories(liblsms PRIVATE 
    ${CMAKE_SOURCE_DIR}/include 
    ${CMAKE_SOURCE_DIR}/src
    ${LUA_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
)

#target_link_libraries(liblsms
#    PUBLIC
#)

set_target_properties(liblsms PROPERTIES OUTPUT_NAME "lsms")

# Targets ------------------------------------------------------------------
add_executable(lsms lsms.cpp)
target_link_libraries(lsms liblsms ${LIBS})

add_executable(lsms_inter EXCLUDE_FROM_ALL lsms_inter.cpp)
target_link_libraries(lsms_inter liblsms)

# needs mjson
add_executable(wl-lsms wl_lsms.cpp)
target_link_libraries(wl-lsms liblsms mjson ${LIBS})

add_executable(rewl-lsms EXCLUDE_FROM_ALL rewl_lsms.cpp)
target_link_libraries(rewl-lsms liblsms)

add_executable(zblock_lu_driver EXCLUDE_FROM_ALL zblock_lu_cuda_driver.cpp)
target_link_libraries(zblock_lu_driver liblsms)
set_target_properties (zblock_lu_driver PROPERTIES COMPILE_DEFINITIONS "NO_HDF5")

# include Directories ----------------------------------------------------------
target_include_directories(lsms PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lua/include
    ${CMAKE_SOURCE_DIR}/src
)

target_include_directories(wl-lsms PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lua/include
    ${CMAKE_SOURCE_DIR}/src
)

target_include_directories(lsms_inter PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lua/include
    ${CMAKE_SOURCE_DIR}/src
)

target_include_directories(rewl-lsms PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lua/include
    ${CMAKE_SOURCE_DIR}/src
)
