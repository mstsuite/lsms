target_sources(
  lsmscore
  PUBLIC calsig.f
         caltnode.f
         celbnd.f
         chkbnd.f
         chkedge.f
         filter_edge.f
         inter_dip.f
         inter.f
         interstitial.f
         intphi.f
         intpl0.f
         invm3.f
         polyhedron.f
         rcritpts.f90
         setup_boundary.f
         setup_vorpol.f
         sigma.f
         sort.f
         sortidx.f
         stepyll.f
         volvor.f
         VORPOL.hpp
         setupVorpol.cpp
         setupVorpol.hpp)

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
