cmake_minimum_required(VERSION 3.0)

set(LIBPDFTL_HEADERS
  TlAtom.h                     
  TlColVectorMatrix.h          
  TlFile.h                     
  TlFileMatrix.h               
  TlFileSymmetricMatrix.h      
  TlFmt.h                      
  TlGetopt.h                   
  TlIntVector.h                
  TlLebedevGrid.h              
  TlListMatrix.h               
  TlLogging.h                  
  TlMath.h                     
  TlMatrix.h                   
  TlMatrixObject.h             
  TlMemManager.h               
  TlMmapAllocater.h            
  TlMsgPack.h                  
  TlPartialMatrix.h            
  TlPartialSymmetricMatrix.h   
  TlPosition.h                 
  TlPrdctbl.h                  
  TlPseudoYaml.h               
  TlResidue.h                  
  TlRowVectorMatrix.h          
  TlSerializeData.h            
  TlSharedPointer.h            
  TlSimpleVector.h             
  TlSparseHashMatrix.h         
  TlSparseHashSymmetricMatrix.h
  TlSparseMatrix.h         
  TlSparseSymmetricMatrix.h
  TlSparseVector.h
  TlSparseVectorMatrix.h
  TlStlUtils.h          
  TlStringTokenizer.h   
  TlSymmetricMatrix.h   
  TlSystem.h            
  TlTime.h              
  TlUtils.h             
  TlVector.h            
  TlVectorMatrixObject.h
  TlVectorObject.h
  )


set(LIBPDFTL_SOURCES 
  TlAtom.cpp
  TlColVectorMatrix.cpp
  TlFile.cpp
  TlFileMatrix.cpp
  TlFileSymmetricMatrix.cpp
  TlFmt.cpp
  TlFname.cxx
  TlGetopt.cpp
  TlLebedevGrid.cpp
  TlListMatrix.cpp
  TlLogging.cpp
  TlMath.cxx
  TlMatrix.cpp
  TlMatrixObject.cpp
  TlMemManager.cpp
  TlMsgPack.cpp
  TlPartialMatrix.cpp
  TlPartialSymmetricMatrix.cpp
  TlPosition.cpp
  TlPrdctbl.cxx
  TlPseudoYaml.cpp
  TlResidue.cxx
  TlRowVectorMatrix.cpp
  TlSerializeData.cpp
  TlSparseHashMatrix.cpp
  TlSparseHashSymmetricMatrix.cpp
  TlSparseMatrix.cpp
  TlSparseSymmetricMatrix.cpp
  TlSparseVector.cpp
  TlSparseVectorMatrix.cpp
  TlStringTokenizer.cpp
  TlSystem.cpp
  TlSymmetricMatrix.cpp  
  TlTime.cpp
  TlUtils.cpp
  TlVector.cpp
  TlVectorMatrixObject.cpp
  )


set(LIBPPDFTL_HEADERS
  TlCommunicate.h
  TlDistributeMatrix.h
  TlDistributeSymmetricMatrix.h
  TlDistributeVector.h
  scalapack.h
  )


set(LIBPPDFTL_SOURCES
  TlCommunicate.cpp
  TlDistributeMatrix.cpp
  TlDistributeSymmetricMatrix.cpp
  TlDistributeVector.cpp
  )


add_library(libpdftl ${LIBPDFTL_HEADERS} ${LIBPDFTL_SOURCES})
target_link_libraries(libpdftl ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
set_target_properties(libpdftl PROPERTIES OUTPUT_NAME "pdftl")

install(TARGETS libpdftl
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
install(FILES ${LIBPDFTL_HEADERS} DESTINATION include)


if (MPI_FOUND)
  add_library(libppdftl ${LIBPPDFTL_HEADERS} ${LIBPPDFTL_SOURCES})
  target_link_libraries(libppdftl ${SCALAPACK_LIBRARIES})
  set_target_properties(libppdftl PROPERTIES OUTPUT_NAME "ppdftl")

  install(TARGETS libppdftl
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
  install(FILES ${LIBPPDFTL_HEADERS} DESTINATION include)
endif(MPI_FOUND)
