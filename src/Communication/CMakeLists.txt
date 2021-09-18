#
# Created by Franco P. Moitzi
#

target_sources(
  lsmscore
  PUBLIC distributeAtoms.cpp distributeAtoms.hpp LSMSCommunication.cpp
         LSMSCommunication.hpp REWLCommunication.cpp REWLCommunication.hpp)

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
