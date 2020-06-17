set(LIBPDFTL_HEADERS
  tl_assert.h
  tl_dense_general_matrix_object.h
  tl_dense_symmetric_matrix_object.h
  tl_dense_matrix_impl_object.h
  tl_dense_vector_object.h
  tl_dense_vector_impl_object.h
  tl_matrix_object.h
  tl_matrix_utils.h
  tl_sparse_general_matrix_object.h
  tl_sparse_symmetric_matrix_object.h
  tl_sparse_matrix_impl_object.h
# ---------------------------------------------------------------------------
  TlAtom.h
  TlFile.h
  TlFmt.h
  TlGetopt.h
  TlIntVector.h
  TlLebedevGrid.h
  TlLogging.h
  TlMath.h
  TlMemManager.h
  TlMmapAllocater.h

  tl_dense_general_matrix_arrays_coloriented.h
  tl_dense_general_matrix_arrays_roworiented.h
  tl_dense_general_matrix_io.h
  tl_dense_general_matrix_mmap.h

  tl_dense_matrix_arrays_object.h
  tl_dense_matrix_io_object.h
  tl_dense_matrix_mmap_object.h

  tl_dense_symmetric_matrix_abstract.h
  tl_dense_symmetric_matrix_io.h

  tl_matrix_mmap_row_major.h
  tl_symmetric_matrix_object.h

  tl_sparse_matrix.h
  tl_sparse_symmetric_matrix.h
  tl_sparse_vector.h

  TlMsgPack.h
  TlPosition.h
  TlPrdctbl.h
  TlPseudoYaml.h
  TlResidue.h
  TlSerializeData.h
  TlSharedPointer.h
  TlStlUtils.h
  TlStringTokenizer.h
  TlSystem.h
  TlTime.h
  TlUtils.h

  tl_vector_object.h
  tl_vector_abstract.h
  tl_vector_utils.h

  # ---------------------------------------------------------------------------
  )


set(LIBPDFTL_SOURCES
  TlAtom.cpp
  tl_dense_general_matrix_arrays_coloriented.cc
  TlFile.cpp
  TlFmt.cpp
  TlFname.cpp
  TlGetopt.cpp
  TlLebedevGrid.cpp
  TlLogging.cpp
  TlMath.cpp
  tl_matrix_object.cc
  tl_matrix_utils.cc
  TlMemManager.cpp

  tl_dense_general_matrix_arrays_roworiented.cc
  tl_dense_general_matrix_io.cc
  tl_dense_general_matrix_mmap.cc
  tl_dense_symmetric_matrix_abstract.cc
  tl_dense_symmetric_matrix_io.cc
  tl_matrix_mmap_row_major.cc

  tl_dense_matrix_arrays_object.cc
  tl_dense_matrix_mmap_object.cc
  tl_dense_matrix_io_object.cc

  tl_sparse_matrix.cc
  tl_sparse_symmetric_matrix.cc
  tl_sparse_vector.cc

  TlMsgPack.cpp
  TlPosition.cpp
  TlPrdctbl.cpp
  TlPseudoYaml.cpp
  TlResidue.cpp
  TlSerializeData.cpp
  TlStringTokenizer.cpp
  TlSystem.cpp
  TlTime.cpp
  TlUtils.cpp

  tl_vector_utils.cc

  # ---------------------------------------------------------------------------
  tl_dense_general_matrix_object.cc
  tl_dense_symmetric_matrix_object.cc
  tl_dense_matrix_impl_object.cc

  tl_dense_vector_object.cc
  tl_dense_vector_impl_object.cc

  tl_sparse_general_matrix_object.cc
  tl_sparse_symmetric_matrix_object.cc
  tl_sparse_matrix_impl_object.cc
  )

# objects using Lapack
if (LAPACK_FOUND)
set(LIBPDFTL_HEADERS
  ${LIBPDFTL_HEADERS}
  lapack.h
  tl_dense_general_matrix_lapack.h
  tl_dense_general_matrix_impl_lapack.h
  tl_dense_symmetric_matrix_lapack.h
  tl_dense_symmetric_matrix_impl_lapack.h
  tl_dense_vector_lapack.h
  tl_dense_vector_impl_lapack.h
  )

set(LIBPDFTL_SOURCES
  ${LIBPDFTL_SOURCES}
  tl_dense_general_matrix_lapack.cc
  tl_dense_general_matrix_impl_lapack.cc
  tl_dense_symmetric_matrix_lapack.cc
  tl_dense_symmetric_matrix_impl_lapack.cc
  tl_dense_vector_lapack.cc
  tl_dense_vector_impl_lapack.cc
  )

set(LIBPDFTL_LIBRARIES
  ${LIBPDFTL_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
endif(LAPACK_FOUND)

# objects using Eigen3
if (EIGEN3_FOUND)
  set(LIBPDFTL_HEADERS
    ${LIBPDFTL_HEADERS}
    tl_dense_general_matrix_eigen.h
    tl_dense_general_matrix_impl_eigen.h
    tl_dense_symmetric_matrix_eigen.h
    tl_dense_symmetric_matrix_impl_eigen.h
    tl_dense_vector_eigen.h
    tl_dense_vector_impl_eigen.h

    tl_sparse_general_matrix_eigen.h
    tl_sparse_general_matrix_impl_eigen.h
    tl_sparse_symmetric_matrix_eigen.h
    tl_sparse_symmetric_matrix_impl_eigen.h

    tl_system_eigen.h
    )

  set(LIBPDFTL_SOURCES
    ${LIBPDFTL_SOURCES}
    tl_dense_general_matrix_eigen.cc
    tl_dense_general_matrix_impl_eigen.cc
    tl_dense_symmetric_matrix_eigen.cc
    tl_dense_symmetric_matrix_impl_eigen.cc
    tl_dense_vector_eigen.cc
    tl_dense_vector_impl_eigen.cc

    tl_sparse_general_matrix_eigen.cc
    tl_sparse_general_matrix_impl_eigen.cc
    tl_sparse_symmetric_matrix_eigen.cc
    tl_sparse_symmetric_matrix_impl_eigen.cc
    )
endif(EIGEN3_FOUND)

# objects using ViennaCL
if (VIENNACL_FOUND)
  set(LIBPDFTL_HEADERS
    ${LIBPDFTL_HEADERS}
    tl_viennacl.h
    tl_dense_general_matrix_viennacl.h
    tl_dense_general_matrix_impl_viennacl.h
    tl_dense_symmetric_matrix_viennacl.h
    tl_dense_symmetric_matrix_impl_viennacl.h
    tl_dense_vector_viennacl.h
    tl_dense_vector_impl_viennacl.h

    tl_sparse_general_matrix_viennacl.h
    tl_sparse_general_matrix_impl_viennacl.h
    tl_sparse_symmetric_matrix_viennacl.h
    tl_sparse_symmetric_matrix_impl_viennacl.h
    )

  set(LIBPDFTL_SOURCES
    ${LIBPDFTL_SOURCES}
    tl_viennacl.cc
    tl_dense_general_matrix_viennacl.cc
    tl_dense_general_matrix_impl_viennacl.cc
    tl_dense_symmetric_matrix_viennacl.cc
    tl_dense_symmetric_matrix_impl_viennacl.cc
    tl_dense_vector_viennacl.cc
    tl_dense_vector_impl_viennacl.cc

    tl_sparse_general_matrix_viennacl.cc
    tl_sparse_general_matrix_impl_viennacl.cc
    tl_sparse_symmetric_matrix_viennacl.cc
    tl_sparse_symmetric_matrix_impl_viennacl.cc
    )

  set(LIBPDFTL_LIBRARIES
    ${LIBPDFTL_LIBRARIES}
    ${OpenCL_LIBRARIES}
  )
endif(VIENNACL_FOUND)

# objects using HDF5
if (HDF5_FOUND)
  set(LIBPDFTL_HEADERS
    ${LIBPDFTL_HEADERS}
    TlHdf5Utils.h
    )

  set(LIBPDFTL_SOURCES
    ${LIBPDFTL_SOURCES}
    TlHdf5Utils.cpp
    )

  set(LIBPDFTL_LIBRARIES
    ${LIBPDFTL_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_LIBRARIES}
    )
endif(HDF5_FOUND)


set(LIBPPDFTL_HEADERS
  TlCommunicate.h
  tl_scalapack_context.h
  tl_dense_scalapack_object.h

  tl_dense_general_matrix_scalapack.h
  tl_dense_symmetric_matrix_scalapack.h
  tl_dense_vector_scalapack.h

  tl_dense_general_matrix_impl_scalapack.h
  tl_dense_symmetric_matrix_impl_scalapack.h
  tl_dense_vector_impl_scalapack.h
  scalapack.h
  )


set(LIBPPDFTL_SOURCES
  TlCommunicate.cpp
  tl_scalapack_context.cc
  tl_dense_scalapack_object.cc

  tl_dense_general_matrix_scalapack.cc
  tl_dense_symmetric_matrix_scalapack.cc
  tl_dense_vector_scalapack.cc

  tl_dense_general_matrix_impl_scalapack.cc
  tl_dense_symmetric_matrix_impl_scalapack.cc
  tl_dense_vector_impl_scalapack.cc
  )


add_library(libpdftl
  ${LIBPDFTL_HEADERS}
  ${LIBPDFTL_SOURCES}
  )
target_link_libraries(libpdftl INTERFACE
  ${LIBPDFTL_LIBRARIES}
)
set_target_properties(libpdftl PROPERTIES OUTPUT_NAME "pdftl")

# include path
if (EIGEN3_FOUND)
target_include_directories(libpdftl PUBLIC ${EIGEN3_INCLUDE_DIR})
endif (EIGEN3_FOUND)

if (VIENNACL_FOUND)
target_include_directories(libpdftl PUBLIC ${VIENNACL_INCLUDE_DIR})
endif (VIENNACL_FOUND)


install(TARGETS libpdftl
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${LIBPDFTL_HEADERS} DESTINATION include)


if (MPI_FOUND)
  add_library(libppdftl
    ${LIBPPDFTL_HEADERS}
    ${LIBPPDFTL_SOURCES}
    )

  if (SCALAPACK_FOUND)
  target_link_libraries(libppdftl INTERFACE
    ${SCALAPACK_LIBRARIES}
    ${BLACS_LIBRARIES}
    ${BLACS_CINIT_LIBRARIES}
    ${BLACS_F77INIT_LIBRARIES}
    )
  endif()
  target_link_libraries(libppdftl INTERFACE
    ${MPI_CXX_LIBRARIES}
    )
  set_target_properties(libppdftl PROPERTIES OUTPUT_NAME "ppdftl")

  install(TARGETS libppdftl
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
  install(FILES ${LIBPPDFTL_HEADERS} DESTINATION include)
endif(MPI_FOUND)
