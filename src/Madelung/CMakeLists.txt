target_sources(
  lsmscore
  PUBLIC cal_madelung_matrix.f
         getkncut.f
         getrscut.f
         getstruc.f
         interf.f
         interfsmr.f
         lattice.f
         madewd.f
         madewdj.f
         madsum.f
         ord3v.f
         pqintg_c.f
         lmfacts.f
         bessj.f
         madelung.h
         Madelung.cpp
         Madelung.hpp)

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
