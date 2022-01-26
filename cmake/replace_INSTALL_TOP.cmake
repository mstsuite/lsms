 message("Input File: " ${CMAKE_ARGV3})
 message("Output File: " ${CMAKE_ARGV4})
 message("New INSTALL_TOP: " ${CMAKE_ARGV5})

file(STRINGS ${CMAKE_ARGV3} inFile)

file(WRITE ${CMAKE_ARGV4})
foreach(line IN LISTS inFile)
  if(${line} MATCHES "^INSTALL_TOP")
    file(APPEND ${CMAKE_ARGV4} "INSTALL_TOP=${CMAKE_ARGV5}\n")
  else()
    file(APPEND ${CMAKE_ARGV4} "${line}\n")
  endif()
endforeach()

