cmake_minimum_required(VERSION 3.5)

set(PDF_COMPONENT_SOURCES
  main_component.cpp
  )
add_executable(pdf-component ${PDF_COMPONENT_SOURCES})
target_link_libraries(pdf-component libpdf libpdftl)
install(TARGETS pdf-component RUNTIME DESTINATION bin)


set(PDF_ERI_SOURCES
  main_eri.cpp
  )
add_executable(pdf-eri ${PDF_ERI_SOURCES})
target_link_libraries(pdf-eri libpdf libpdftl)
install(TARGETS pdf-eri RUNTIME DESTINATION bin)


set(PDF_I2PQ_SHOW_SOURCES
  main_i2pq_show.cpp
  )
add_executable(pdf-i2pq-show ${PDF_I2PQ_SHOW_SOURCES})
target_link_libraries(pdf-i2pq-show libpdf libpdftl)
install(TARGETS pdf-i2pq-show RUNTIME DESTINATION bin)


set(PDF_INFO_BASISSET_SOURCES
  main_info_basisset.cpp
  )
add_executable(pdf-info-basisset ${PDF_INFO_BASISSET_SOURCES})
target_link_libraries(pdf-info-basisset libpdf libpdftl)
install(TARGETS pdf-info-basisset RUNTIME DESTINATION bin)


set(PDF_INFO_CHARGE_SOURCES
  main_info_charge.cpp
  )
add_executable(pdf-info-charge ${PDF_INFO_CHARGE_SOURCES})
target_link_libraries(pdf-info-charge libpdf libpdftl)
install(TARGETS pdf-info-charge RUNTIME DESTINATION bin)


set(PDF_INFO_GRID_SOURCES
  main_info_grid.cpp
  )
add_executable(pdf-info-grid ${PDF_INFO_GRID_SOURCES})
target_link_libraries(pdf-info-grid libpdf libpdftl)
install(TARGETS pdf-info-grid RUNTIME DESTINATION bin)


set(PDF_INFO_ORB_SOURCES
  main_info_orb.cpp
  )
add_executable(pdf-info-orb ${PDF_INFO_ORB_SOURCES})
target_link_libraries(pdf-info-orb libpdf libpdftl)
install(TARGETS pdf-info-orb RUNTIME DESTINATION bin)


set(PDF_INIT_PARAM_SOURCES
  main_init_param.cpp
  )
add_executable(pdf-init-param ${PDF_INIT_PARAM_SOURCES})
target_link_libraries(pdf-init-param libpdf libpdftl)
install(TARGETS pdf-init-param RUNTIME DESTINATION bin)


set(PDF_KWD_SOURCES
  main_kwd.cpp
  )
add_executable(pdf-kwd ${PDF_KWD_SOURCES})
target_link_libraries(pdf-kwd libpdf libpdftl)
install(TARGETS pdf-kwd RUNTIME DESTINATION bin)


set(PDF_LO_SOURCES
  main_lo.cpp
  DfLocalize.h
  DfLocalize.cpp
  )
add_executable(pdf-lo ${PDF_LO_SOURCES})
target_link_libraries(pdf-lo libpdf libpdftl)
install(TARGETS pdf-lo RUNTIME DESTINATION bin)


set(PDF_MAT_ADD_SOURCES
  main_mat_add.cpp
  )
add_executable(pdf-mat-add ${PDF_MAT_ADD_SOURCES})
target_link_libraries(pdf-mat-add libpdf libpdftl)
install(TARGETS pdf-mat-add RUNTIME DESTINATION bin)


set(PDF_MAT_DIAGONAL_SOURCES
  main_mat_diagonal.cpp
  )
add_executable(pdf-mat-diagonal ${PDF_MAT_DIAGONAL_SOURCES})
target_link_libraries(pdf-mat-diagonal libpdf libpdftl)
install(TARGETS pdf-mat-diagonal RUNTIME DESTINATION bin)


set(PDF_MAT_DIFF_SOURCES
  main_mat_diff.cpp
  )
add_executable(pdf-mat-diff ${PDF_MAT_DIFF_SOURCES})
target_link_libraries(pdf-mat-diff libpdf libpdftl)
install(TARGETS pdf-mat-diff RUNTIME DESTINATION bin)


set(PDF_MAT_DOT_SOURCES
  main_mat_dot.cpp
  )
add_executable(pdf-mat-dot ${PDF_MAT_DOT_SOURCES})
target_link_libraries(pdf-mat-dot libpdf libpdftl)
install(TARGETS pdf-mat-dot RUNTIME DESTINATION bin)


set(PDF_MAT_EXTENT_SOURCES
  main_mat_extent.cpp
  )
add_executable(pdf-mat-extent ${PDF_MAT_EXTENT_SOURCES})
target_link_libraries(pdf-mat-extent libpdf libpdftl)
install(TARGETS pdf-mat-extent RUNTIME DESTINATION bin)


set(PDF_MAT_INFO_SOURCES
  main_mat_info.cpp
  )
add_executable(pdf-mat-info ${PDF_MAT_INFO_SOURCES})
target_link_libraries(pdf-mat-info libpdf libpdftl)
install(TARGETS pdf-mat-info RUNTIME DESTINATION bin)


set(PDF_MAT_MUL_SOURCES
  main_mat_mul.cpp
  )
add_executable(pdf-mat-mul ${PDF_MAT_MUL_SOURCES})
target_link_libraries(pdf-mat-mul libpdf libpdftl)
install(TARGETS pdf-mat-mul RUNTIME DESTINATION bin)


set(PDF_MAT_RESIZE_SOURCES
  main_mat_resize.cpp
  )
add_executable(pdf-mat-resize ${PDF_MAT_RESIZE_SOURCES})
target_link_libraries(pdf-mat-resize libpdf libpdftl)
install(TARGETS pdf-mat-resize RUNTIME DESTINATION bin)


set(PDF_MAT_SELECT_SOURCES
  main_mat_select.cpp
  )
add_executable(pdf-mat-select ${PDF_MAT_SELECT_SOURCES})
target_link_libraries(pdf-mat-select libpdf libpdftl)
install(TARGETS pdf-mat-select RUNTIME DESTINATION bin)


set(PDF_MAT_SHOW_SOURCES
  main_mat_show.cpp
  )
add_executable(pdf-mat-show ${PDF_MAT_SHOW_SOURCES})
target_link_libraries(pdf-mat-show libpdf libpdftl)
install(TARGETS pdf-mat-show RUNTIME DESTINATION bin)


set(PDF_MAT_SYMMETRIZE_SOURCES
  main_mat_symmetrize.cpp
  )
add_executable(pdf-mat-symmetrize ${PDF_MAT_SYMMETRIZE_SOURCES})
target_link_libraries(pdf-mat-symmetrize libpdf libpdftl)
install(TARGETS pdf-mat-symmetrize RUNTIME DESTINATION bin)


set(PDF_MAT_TRANSPOSE_SOURCES
  main_mat_transpose.cpp
  )
add_executable(pdf-mat-transpose ${PDF_MAT_TRANSPOSE_SOURCES})
target_link_libraries(pdf-mat-transpose libpdf libpdftl)
install(TARGETS pdf-mat-transpose RUNTIME DESTINATION bin)


set(PDF_MKFLD_MO_SOURCES
  main_mkfld_mo.cpp
  mkfld_common.h
  mkfld_common.cpp
  TlFieldDataObject.h
  TlFieldDataObject.cpp
  TlMoField.h
  TlMoField.cpp
  )
add_executable(pdf-mkfld-mo ${PDF_MKFLD_MO_SOURCES})
target_link_libraries(pdf-mkfld-mo libpdf libpdftl)
install(TARGETS pdf-mkfld-mo RUNTIME DESTINATION bin)


set(PDF_MKFLD_DENS_SOURCES
  main_mkfld_dens.cpp
  mkfld_common.h
  mkfld_common.cpp
  TlFieldDataObject.h
  TlFieldDataObject.cpp
  TlDensField.h
  TlDensField.cpp
  )
add_executable(pdf-mkfld-dens ${PDF_MKFLD_DENS_SOURCES})
target_link_libraries(pdf-mkfld-dens libpdf libpdftl)
install(TARGETS pdf-mkfld-dens RUNTIME DESTINATION bin)


set(PDF_MKFLD_ESP_SOURCES
  main_mkfld_esp.cpp
  mkfld_common.h
  mkfld_common.cpp
  TlFieldDataObject.h
  TlFieldDataObject.cpp
  TlEspField.h
  TlEspField.cpp
  )
add_executable(pdf-mkfld-esp ${PDF_MKFLD_ESP_SOURCES})
target_link_libraries(pdf-mkfld-esp libpdf libpdftl)
install(TARGETS pdf-mkfld-esp RUNTIME DESTINATION bin)


set(PDF_POP_MULLIKEN_SOURCES
  main_pop_mulliken.cpp
  )
add_executable(pdf-pop-mulliken ${PDF_POP_MULLIKEN_SOURCES})
target_link_libraries(pdf-pop-mulliken libpdf libpdftl)
install(TARGETS pdf-pop-mulliken RUNTIME DESTINATION bin)


set(PDF_POP_ESP_SOURCES
  main_pop_esp.cpp
  TlFieldDataObject.h
  TlFieldDataObject.cpp
  TlEspField.h
  TlEspField.cpp
  )
add_executable(pdf-pop-esp ${PDF_POP_ESP_SOURCES})
target_link_libraries(pdf-pop-esp libpdf libpdftl)
install(TARGETS pdf-pop-esp RUNTIME DESTINATION bin)


set(PDF_VTR_DIFF_SOURCES
  main_vtr_diff.cpp
  )
add_executable(pdf-vtr-diff ${PDF_VTR_DIFF_SOURCES})
target_link_libraries(pdf-vtr-diff libpdf libpdftl)
install(TARGETS pdf-vtr-diff RUNTIME DESTINATION bin)


set(PDF_VTR_SHOW_SOURCES
  main_vtr_show.cpp
  )
add_executable(pdf-vtr-show ${PDF_VTR_DIFF_SOURCES})
target_link_libraries(pdf-vtr-show libpdf libpdftl)
install(TARGETS pdf-vtr-show RUNTIME DESTINATION bin)


if (MPI_FOUND)
  set(PDF_LO_MPI_SOURCES
    main_lo_mpi.cpp
    DfLocalize.h
    DfLocalize.cpp
    DfLocalize_Parallel.h
    DfLocalize_Parallel.cpp
    )
  add_executable(pdf-lo-mpi ${PDF_LO_MPI_SOURCES})
  target_link_libraries(pdf-lo-mpi libppdf libppdftl libpdf libpdftl ${MPI_CXX_LIBRARIES})
  install(TARGETS pdf-lo-mpi RUNTIME DESTINATION bin)
endif(MPI_FOUND)
