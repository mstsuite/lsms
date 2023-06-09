

target_sources(lsmscore
        PUBLIC

        ChargeDensity.cpp
        ChargeDensity.hpp

        )

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)