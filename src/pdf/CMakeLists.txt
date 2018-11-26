# -----------------------------------------------------------------------------
# for PDF.x
# -----------------------------------------------------------------------------
set(PDF_SOURCES
  main.cpp
  )

add_executable(PDF.x ${PDF_SOURCES})
target_link_libraries(PDF.x
  libpdf
  libpdftl
  )
install(TARGETS PDF.x RUNTIME DESTINATION bin)


# -----------------------------------------------------------------------------
# for PPDF.x
# -----------------------------------------------------------------------------
set(PPDF_SOURCES
  main_Parallel.cpp
  )

if (MPI_FOUND)
  add_executable(PPDF.x ${PPDF_SOURCES})
  target_link_libraries(PPDF.x
    libppdf
    libpdf
    libppdftl
    libpdftl
    )
  install(TARGETS PPDF.x RUNTIME DESTINATION bin)
endif()
