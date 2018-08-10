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
  ${OpenCL_LIBRARIES}
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
    libppdftl 
    libpdf 
    libpdftl
    ${HDF5_CXX_LIBRARIES}  
    ${MPI_CXX_LIBRARIES} 
    ${OpenCL_LIBRARIES})
  install(TARGETS PPDF.x RUNTIME DESTINATION bin)
endif()
