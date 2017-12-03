cmake_minimum_required(VERSION 3.5)

set(PDF_SOURCES
  main.cpp
  )

set(PPDF_SOURCES
  main_Parallel.cpp
  )

add_executable(PDF.x ${PDF_SOURCES})
target_link_libraries(PDF.x libpdf libpdftl)
install(TARGETS PDF.x RUNTIME DESTINATION bin)

if (MPI_FOUND)
  add_executable(PPDF.x ${PPDF_SOURCES})
  target_link_libraries(PPDF.x libppdf libppdftl libpdf libpdftl ${MPI_CXX_LIBRARIES})
  install(TARGETS PPDF.x RUNTIME DESTINATION bin)
endif()
