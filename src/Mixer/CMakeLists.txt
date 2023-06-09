target_sources(lsmscore
        PUBLIC

        MixingParameter.hpp
        MixingParameter.cpp

        Mixer.hpp
        Mixer.cpp

        MixerType.hpp
        MixerType.cpp

        MixingVector.hpp
        MixingVector.cpp
        )

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)