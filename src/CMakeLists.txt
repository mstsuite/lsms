
target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(Accelerator)

add_subdirectory(Communication)

add_subdirectory(Core)

add_subdirectory(Kubo)

add_subdirectory(ChargeDensity)

add_subdirectory(LuaInterface)

add_subdirectory(Main)

add_subdirectory(Madelung)

add_subdirectory(MultipoleMadelung)

add_subdirectory(Mixer)

add_subdirectory(Misc)

add_subdirectory(MultipleScattering)

add_subdirectory(Potential)

add_subdirectory(RadialGrid)

add_subdirectory(SingleSite)

add_subdirectory(TotalEnergy)

add_subdirectory(VORPOL)

