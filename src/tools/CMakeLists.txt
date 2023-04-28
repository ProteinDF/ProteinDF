set(PDF_COMPONENT_SOURCES main_component.cpp)
add_executable(pdf-component ${PDF_COMPONENT_SOURCES})
target_link_libraries(pdf-component libpdf libpdftl)
install(TARGETS pdf-component RUNTIME DESTINATION bin)

set(PDF_ERI_SOURCES main_eri.cpp)
add_executable(pdf-eri ${PDF_ERI_SOURCES})
target_link_libraries(pdf-eri libpdf libpdftl)
install(TARGETS pdf-eri RUNTIME DESTINATION bin)

set(PDF_I2PQ_SHOW_SOURCES main_i2pq_show.cpp)
add_executable(pdf-i2pq-show ${PDF_I2PQ_SHOW_SOURCES})
target_link_libraries(pdf-i2pq-show libpdf libpdftl)
install(TARGETS pdf-i2pq-show RUNTIME DESTINATION bin)

set(PDF_INFO_BASISSET_SOURCES main_info_basisset.cpp)
add_executable(pdf-info-basisset ${PDF_INFO_BASISSET_SOURCES})
target_link_libraries(pdf-info-basisset libpdf libpdftl)
install(TARGETS pdf-info-basisset RUNTIME DESTINATION bin)

set(PDF_INFO_CHARGE_SOURCES main_info_charge.cpp)
add_executable(pdf-info-charge ${PDF_INFO_CHARGE_SOURCES})
target_link_libraries(pdf-info-charge libpdf libpdftl)
install(TARGETS pdf-info-charge RUNTIME DESTINATION bin)

set(PDF_INFO_MO_OVERLAP_SOURCES main_info_mo_overlap.cc)
add_executable(pdf-info-mo-overlap ${PDF_INFO_MO_OVERLAP_SOURCES})
target_link_libraries(pdf-info-mo-overlap libpdf libpdftl)
install(TARGETS pdf-info-mo-overlap RUNTIME DESTINATION bin)

set(PDF_INFO_ECOMP_SOURCES main_info_ecomponents.cc)
add_executable(pdf-info-ecomp ${PDF_INFO_ECOMP_SOURCES})
target_link_libraries(pdf-info-ecomp libpdf libpdftl)
install(TARGETS pdf-info-ecomp RUNTIME DESTINATION bin)

set(PDF_INFO_GRID_SOURCES main_info_grid.cpp)
add_executable(pdf-info-grid ${PDF_INFO_GRID_SOURCES})
target_link_libraries(pdf-info-grid libpdf libpdftl)
install(TARGETS pdf-info-grid RUNTIME DESTINATION bin)

set(PDF_INFO_ORB_SOURCES main_info_orb.cpp)
add_executable(pdf-info-orb ${PDF_INFO_ORB_SOURCES})
target_link_libraries(pdf-info-orb libpdf libpdftl)
install(TARGETS pdf-info-orb RUNTIME DESTINATION bin)

if(VIENNACL_FOUND)
  set(PDF_INFO_VIENNACL_SOURCES main_info_viennacl.cc)
  add_executable(pdf-info-viennacl ${PDF_INFO_VIENNACL_SOURCES})
  target_include_directories(pdf-info-viennacl PUBLIC ${VIENNACL_INCLUDE_DIR})
  target_link_libraries(pdf-info-viennacl ${OpenCL_LIBRARIES})
  install(TARGETS pdf-info-viennacl RUNTIME DESTINATION bin)
endif(VIENNACL_FOUND)

set(PDF_INIT_PARAM_SOURCES main_init_param.cpp)
add_executable(pdf-init-param ${PDF_INIT_PARAM_SOURCES})
target_link_libraries(pdf-init-param libpdf libpdftl)
install(TARGETS pdf-init-param RUNTIME DESTINATION bin)

set(PDF_KWD_SOURCES main_kwd.cpp)
add_executable(pdf-kwd ${PDF_KWD_SOURCES})
target_link_libraries(pdf-kwd libpdf libpdftl)
install(TARGETS pdf-kwd RUNTIME DESTINATION bin)

set(PDF_LO_SOURCES main_lo.cpp DfLocalize.h DfLocalize.cpp)
add_executable(pdf-lo ${PDF_LO_SOURCES})
target_link_libraries(pdf-lo libpdf libpdftl)
install(TARGETS pdf-lo RUNTIME DESTINATION bin)

set(PDF_LO_1_SOURCES main_lo.1.cc DfLocalize.1.h DfLocalize.1.cc)
add_executable(pdf-lo.1 ${PDF_LO_1_SOURCES})
target_link_libraries(pdf-lo.1 libpdf libpdftl)
install(TARGETS pdf-lo.1 RUNTIME DESTINATION bin)

# install(
#   CODE "execute_process( \
#     COMMAND ${CMAKE_COMMAND} -E create_symlink \
#     ${CMAKE_INSTALL_PREFIX}/bin/pdf-lo.1 \
#     ${CMAKE_INSTALL_PREFIX}/bin/pdf-lo \
#     )")

set(PDF_MAT_ADD_SOURCES main_mat_add.cpp)
add_executable(pdf-mat-add ${PDF_MAT_ADD_SOURCES})
target_link_libraries(pdf-mat-add libpdf libpdftl)
install(TARGETS pdf-mat-add RUNTIME DESTINATION bin)

set(PDF_MAT_CHOLESKY_SOURCES main_mat_cholesky.cpp)
add_executable(pdf-mat-cholesky ${PDF_MAT_CHOLESKY_SOURCES})
target_link_libraries(pdf-mat-cholesky libpdf libpdftl)
install(TARGETS pdf-mat-cholesky RUNTIME DESTINATION bin)

set(PDF_MAT_DIAGONAL_SOURCES main_mat_diagonal.cpp)
add_executable(pdf-mat-diagonal ${PDF_MAT_DIAGONAL_SOURCES})
target_link_libraries(pdf-mat-diagonal libpdf libpdftl)
install(TARGETS pdf-mat-diagonal RUNTIME DESTINATION bin)

set(PDF_MAT_DIFF_SOURCES main_mat_diff.cpp)
add_executable(pdf-mat-diff ${PDF_MAT_DIFF_SOURCES})
target_link_libraries(pdf-mat-diff libpdf libpdftl)
install(TARGETS pdf-mat-diff RUNTIME DESTINATION bin)

set(PDF_MAT_DOT_SOURCES main_mat_dot.cpp)
add_executable(pdf-mat-dot ${PDF_MAT_DOT_SOURCES})
target_link_libraries(pdf-mat-dot libpdf libpdftl)
install(TARGETS pdf-mat-dot RUNTIME DESTINATION bin)

set(PDF_MAT_EXTEND_SOURCES main_mat_extend.cc)
add_executable(pdf-mat-extend ${PDF_MAT_EXTEND_SOURCES})
target_link_libraries(pdf-mat-extend libpdf libpdftl)
install(TARGETS pdf-mat-extend RUNTIME DESTINATION bin)

set(PDF_MAT_INFO_SOURCES main_mat_info.cpp)
add_executable(pdf-mat-info ${PDF_MAT_INFO_SOURCES})
target_link_libraries(pdf-mat-info libpdf libpdftl)
install(TARGETS pdf-mat-info RUNTIME DESTINATION bin)


set(PDF_MAT_ISFINITE_SOURCES
  main_mat_isfinite.cc
  )
add_executable(pdf-mat-isfinite ${PDF_MAT_ISFINITE_SOURCES})
target_link_libraries(pdf-mat-isfinite libpdf libpdftl)
install(TARGETS pdf-mat-isfinite RUNTIME DESTINATION bin)


set(PDF_MAT_MUL_SOURCES
  main_mat_mul.cpp
  )
add_executable(pdf-mat-mul ${PDF_MAT_MUL_SOURCES})
target_link_libraries(pdf-mat-mul libpdf libpdftl)
install(TARGETS pdf-mat-mul RUNTIME DESTINATION bin)

set(PDF_MAT_RESIZE_SOURCES main_mat_resize.cpp)
add_executable(pdf-mat-resize ${PDF_MAT_RESIZE_SOURCES})
target_link_libraries(pdf-mat-resize libpdf libpdftl)
install(TARGETS pdf-mat-resize RUNTIME DESTINATION bin)

set(PDF_MAT_SELECT_SOURCES main_mat_select.cpp)
add_executable(pdf-mat-select ${PDF_MAT_SELECT_SOURCES})
target_link_libraries(pdf-mat-select libpdf libpdftl)
install(TARGETS pdf-mat-select RUNTIME DESTINATION bin)

set(PDF_MAT_SHOW_SOURCES main_mat_show.cpp)
add_executable(pdf-mat-show ${PDF_MAT_SHOW_SOURCES})
target_link_libraries(pdf-mat-show libpdf libpdftl)
install(TARGETS pdf-mat-show RUNTIME DESTINATION bin)

set(PDF_MAT_SUM_SOURCES main_mat_sum.cpp)
add_executable(pdf-mat-sum ${PDF_MAT_SUM_SOURCES})
target_link_libraries(pdf-mat-sum libpdf libpdftl)
install(TARGETS pdf-mat-sum RUNTIME DESTINATION bin)

set(PDF_MAT_SYMMETRIZE_SOURCES main_mat_symmetrize.cpp)
add_executable(pdf-mat-symmetrize ${PDF_MAT_SYMMETRIZE_SOURCES})
target_link_libraries(pdf-mat-symmetrize libpdf libpdftl)
install(TARGETS pdf-mat-symmetrize RUNTIME DESTINATION bin)

set(PDF_MAT_TOP_SOURCES main_mat_top.cpp)
add_executable(pdf-mat-top ${PDF_MAT_TOP_SOURCES})
target_link_libraries(pdf-mat-top libpdf libpdftl)
install(TARGETS pdf-mat-top RUNTIME DESTINATION bin)

set(PDF_MAT_TRANSPOSE_SOURCES main_mat_transpose.cpp)
add_executable(pdf-mat-transpose ${PDF_MAT_TRANSPOSE_SOURCES})
target_link_libraries(pdf-mat-transpose libpdf libpdftl)
install(TARGETS pdf-mat-transpose RUNTIME DESTINATION bin)

# set(PDF_MAT_UPDATE_RVM_SOURCES main_mat_update_rvm.cc )
# add_executable(pdf-mat-update-rvm ${PDF_MAT_UPDATE_RVM_SOURCES})
# target_link_libraries(pdf-mat-update-rvm libpdf libpdftl) install(TARGETS
# pdf-mat-update-rvm RUNTIME DESTINATION bin)

set(PDF_MAT_UPGRADE_RVM_SOURCES main_mat_upgrade_rvm.cc)
add_executable(pdf-mat-upgrade-rvm ${PDF_MAT_UPGRADE_RVM_SOURCES})
target_link_libraries(pdf-mat-upgrade-rvm libpdf libpdftl)
install(TARGETS pdf-mat-upgrade-rvm RUNTIME DESTINATION bin)

set(PDF_MINT_SOURCES main_mint.cpp tl_mint_object.h tl_mint_object.cpp
                     tl_mint_pdf.h tl_mint_pdf.cpp)
set(PDF_MINT_LIBRARIES libpdf libpdftl)
if(LIBCINT_FOUND)
  set(PDF_MINT_SOURCES ${PDF_MINT_SOURCES} tl_mint_libcint.h
                       tl_mint_libcint.cpp)
  set(PDF_MINT_LIBRARIES ${PDF_MINT_LIBRARIES} ${LIBCINT_LIBRARIES})
endif()
add_executable(pdf-mint ${PDF_MINT_SOURCES})
target_link_libraries(pdf-mint ${PDF_MINT_LIBRARIES})
install(TARGETS pdf-mint RUNTIME DESTINATION bin)

set(PDF_MKFLD_MO_SOURCES
    main_mkfld_mo.cpp
    mkfld_common.h
    mkfld_common.cpp
    TlFieldDataObject.h
    TlFieldDataObject.cpp
    TlMoField.h
    TlMoField.cpp)
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
    TlDensField.cpp)
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
    TlEspField.cpp)
add_executable(pdf-mkfld-esp ${PDF_MKFLD_ESP_SOURCES})
target_link_libraries(pdf-mkfld-esp libpdf libpdftl)
install(TARGETS pdf-mkfld-esp RUNTIME DESTINATION bin)


set(PDF_MKFLD_ESP_CLASSIC_SOURCES
  main_mkfld_esp_classic.cpp
  mkfld_common.h
  mkfld_common.cpp
  TlFieldDataObject.h
  TlFieldDataObject.cpp
  TlEspField.h
  TlEspField.cpp
  )
add_executable(pdf-mkfld-esp-classic ${PDF_MKFLD_ESP_CLASSIC_SOURCES})
target_link_libraries(pdf-mkfld-esp-classic libpdf libpdftl)
install(TARGETS pdf-mkfld-esp-classic RUNTIME DESTINATION bin)


set(PDF_POP_MULLIKEN_SOURCES
  main_pop_mulliken.cpp
  )
add_executable(pdf-pop-mulliken ${PDF_POP_MULLIKEN_SOURCES})
target_link_libraries(pdf-pop-mulliken libpdf libpdftl)
install(TARGETS pdf-pop-mulliken RUNTIME DESTINATION bin)

set(PDF_POP_ESP_SOURCES main_pop_esp.cpp TlEspPop.h TlEspPop.cpp TlEspField.h
                        TlEspField.cpp)
add_executable(pdf-pop-esp ${PDF_POP_ESP_SOURCES})
target_link_libraries(pdf-pop-esp libpdf libpdftl)
install(TARGETS pdf-pop-esp RUNTIME DESTINATION bin)

set(PDF_POP_RESP_SOURCES main_pop_resp.cpp TlEspPop.h TlEspPop.cpp TlEspField.h
                         TlEspField.cpp)
add_executable(pdf-pop-resp ${PDF_POP_RESP_SOURCES})
target_link_libraries(pdf-pop-resp libpdf libpdftl)
install(TARGETS pdf-pop-resp RUNTIME DESTINATION bin)

set(PDF_VTR_DIFF_SOURCES main_vtr_diff.cpp)
add_executable(pdf-vtr-diff ${PDF_VTR_DIFF_SOURCES})
target_link_libraries(pdf-vtr-diff libpdf libpdftl)
install(TARGETS pdf-vtr-diff RUNTIME DESTINATION bin)

set(PDF_VTR_SHOW_SOURCES main_vtr_show.cpp)
add_executable(pdf-vtr-show ${PDF_VTR_SHOW_SOURCES})
target_link_libraries(pdf-vtr-show libpdf libpdftl)
install(TARGETS pdf-vtr-show RUNTIME DESTINATION bin)

# ================
set(PDF_CVM2RSFD_SOURCES main_cvm2rsfd.cpp)
add_executable(pdf-cvm2rsfd ${PDF_CVM2RSFD_SOURCES})
target_link_libraries(pdf-cvm2rsfd libpdf libpdftl)
install(TARGETS pdf-cvm2rsfd RUNTIME DESTINATION bin)

set(PDF_CVM2CSFD_SOURCES main_cvm2csfd.cpp)
add_executable(pdf-cvm2csfd ${PDF_CVM2CSFD_SOURCES})
target_link_libraries(pdf-cvm2csfd libpdf libpdftl)
install(TARGETS pdf-cvm2csfd RUNTIME DESTINATION bin)

set(PDF_RVM2CSFD_SOURCES main_rvm2csfd.cpp)
add_executable(pdf-rvm2csfd ${PDF_RVM2CSFD_SOURCES})
target_link_libraries(pdf-rvm2csfd libpdf libpdftl)
install(TARGETS pdf-rvm2csfd RUNTIME DESTINATION bin)

# ------------------------------------------------------------------------------
# benchmark
# ------------------------------------------------------------------------------
set(PDF_BENCH_MAT_MUL_SOURCES main_bench_mat_common.h main_bench_mat_mul.cc)
add_executable(pdf-bench-mat-mul ${PDF_BENCH_MAT_MUL_SOURCES})
target_link_libraries(pdf-bench-mat-mul libpdf libpdftl)
install(TARGETS pdf-bench-mat-mul RUNTIME DESTINATION bin)

set(PDF_BENCH_MAT_DIAG_SOURCES main_bench_mat_common.h main_bench_mat_diag.cc)
add_executable(pdf-bench-mat-diag ${PDF_BENCH_MAT_DIAG_SOURCES})
target_link_libraries(pdf-bench-mat-diag libpdf libpdftl)
install(TARGETS pdf-bench-mat-diag RUNTIME DESTINATION bin)

# ------------------------------------------------------------------------------
# debug
# ------------------------------------------------------------------------------
if(VIENNACL_FOUND)
  set(PDF_DEBUG_MAT_MUL_SOURCES tool_common.h main_debug_mat_mul.cc)
  add_executable(pdf-debug-mat-mul ${PDF_DEBUG_MAT_MUL_SOURCES})
  target_link_libraries(pdf-debug-mat-mul libpdf libpdftl)
  install(TARGETS pdf-debug-mat-mul RUNTIME DESTINATION bin)
endif(VIENNACL_FOUND)

# ================

if(HDF5_FOUND)
  set(PDF_MAT2H5_SOURCES main_mat2h5.cpp)
  add_executable(pdf-mat2h5 ${PDF_MAT2H5_SOURCES})
  target_link_libraries(pdf-mat2h5 libpdf libpdftl)
  install(TARGETS pdf-mat2h5 RUNTIME DESTINATION bin)
endif(HDF5_FOUND)

if(MPI_FOUND)
  set(PDF_LO_1_MPI_SOURCES main_lo.1.mpi.cc
    DfLocalize.1.h DfLocalize.1.cc DfLocalize_Parallel.1.h DfLocalize_Parallel.1.cc)
  add_executable(pdf-lo.1.mpi ${PDF_LO_1_MPI_SOURCES})
  target_link_libraries(pdf-lo.1.mpi libppdf libppdftl libpdf libpdftl)
  install(TARGETS pdf-lo.1.mpi RUNTIME DESTINATION bin)

  # install(
  #   CODE "execute_process( \
  #     COMMAND ${CMAKE_COMMAND} -E create_symlink \
  #     ${CMAKE_INSTALL_PREFIX}/bin/pdf-lo.1.mpi \
  #     ${CMAKE_INSTALL_PREFIX}/bin/pdf-lo.mpi \
  #     )")

  set(PDF_LO_MPI_SOURCES main_lo.mpi.cc
    DfLocalize.h DfLocalize.cpp DfLocalize_Parallel.h DfLocalize_Parallel.cpp)
  add_executable(pdf-lo.mpi ${PDF_LO_MPI_SOURCES})
  target_link_libraries(pdf-lo.mpi libppdf libppdftl libpdf libpdftl)
  install(TARGETS pdf-lo.mpi RUNTIME DESTINATION bin)

  set(PDF_CHECK_MPI_SOURCES main_check.mpi.cc)
  add_executable(pdf-check.mpi ${PDF_CHECK_MPI_SOURCES})
  target_link_libraries(pdf-check.mpi libppdf libppdftl libpdf libpdftl)
  install(TARGETS pdf-check.mpi RUNTIME DESTINATION bin)
endif(MPI_FOUND)
