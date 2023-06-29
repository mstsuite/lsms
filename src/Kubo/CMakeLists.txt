
target_sources(
  lsmscore
  PUBLIC CurrentMatrix.cpp
         Conductivity.cpp
         CurrentMatrix.hpp
         CurrentMatrix.hpp)

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
