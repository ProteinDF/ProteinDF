if PARALLEL
lib_LTLIBRARIES = libpdf.la libppdf.la
else
lib_LTLIBRARIES = libpdf.la
endif

libpdf_la_SOURCES = \
	DfXCFunctional.h		DfXCFunctional.cpp			\
	DfCalcGridX.h			DfCalcGridX.cpp				\
	DfFunctional.h			DfFunctional.cpp			\
	DfFunctional_Slater.h	DfFunctional_Slater.cpp		\
	DfFunctional_Becke88.h	DfFunctional_Becke88.cpp	\
	DfFunctional_Becke88_ExceptLDA.h					\
	DfFunctional_Becke88_ExceptLDA.cpp					\
	DfFunctional_VWN.h		DfFunctional_VWN.cpp		\
	DfFunctional_VWN3.h		DfFunctional_VWN3.cpp		\
	DfFunctional_LYP.h		DfFunctional_LYP.cpp		\
	DfFunctional_HFS.h		DfFunctional_HFS.cpp		\
	DfFunctional_SVWN.h		DfFunctional_SVWN.cpp		\
	DfFunctional_B88LYP.h	DfFunctional_B88LYP.cpp		\
	DfFunctional_B3LYP.h	DfFunctional_B3LYP.cpp		\
	DfFunctional_PW91X.h	DfFunctional_PW91X.cpp		\
	CnError.h				CnError.cxx 				\
	DfObject.h				DfObject.cpp				\
	DfDiffDensityMatrix.h	DfDiffDensityMatrix.cpp		\
	DfCalcGrid.h			DfCalcGrid.cpp				\
	DfCD.h					DfCD.cpp					\
	DfConvcheck.h			DfConvcheck.cxx				\
	DfConverge.h			DfConverge.cpp				\
	DfConverge_Anderson.h	DfConverge_Anderson.cpp		\
	DfConverge_Damping.h	DfConverge_Damping.cpp		\
	DfConverge_DIIS.h		DfConverge_DIIS.cpp			\
	DfCqclomatrix.h			DfCqclomatrix.cxx			\
	DfDensityFittingObject.h							\
	DfDiagonal.h			DfDiagonal.cxx				\
	DfDmatrix.h				DfDmatrix.cxx				\
	DfEngineObject.h									\
	DfJMatrix.h				DfJMatrix.cpp				\
	DfKMatrix.h				DfKMatrix.cpp				\
	DfFockMatrix.h			DfFockMatrix.cpp			\
	DfGenerateGrid.h		DfGenerateGrid.cpp			\
	DfGridFreeXC.h			DfGridFreeXC.cpp			\
	DfInitialGuess.h		DfInitialGuess.cpp			\
	DfInitialGuessHuckel.h	DfInitialGuessHuckel.cpp	\
	DfInitialGuessHarris.h	DfInitialGuessHarris.cpp	\
	DfInputdata.h			DfInputdata.cxx				\
	DfIntegrals.h			DfIntegrals.cxx				\
	DfInvMatrix.h			DfInvMatrix.cpp				\
	DfLevelshift.h			DfLevelshift.cxx			\
	DfPreScf.h				DfPreScf.cpp				\
	DfPopulation.h			DfPopulation.cxx			\
	DfQclo.h				DfQclo.cxx					\
	DfScf.h					DfScf.cxx					\
	DfDensityFittingX.h		DfDensityFittingX.cpp		\
	DfSummary.h				DfSummary.cxx				\
	DfThreeindexintegrals.h	DfThreeindexintegrals.cxx	\
	DfTotalEnergy.h			DfTotalEnergy.cpp			\
	DfTransFmatrix.h		DfTransFmatrix.cxx			\
	DfTransatob.h			DfTransatob.cxx				\
	DfXcenefitting.h		DfXcenefitting.cxx			\
	DfXMatrix.h				DfXMatrix.cpp				\
	DfCleanup.h				DfCleanup.cpp				\
	Fl_Db_Basis.h			Fl_Db_Basis.cxx				\
	Fl_Fragment.h			Fl_Fragment.cxx				\
	Fl_Geometry.h			Fl_Geometry.cxx				\
	Fl_Gto.h				Fl_Gto.cxx					\
	Fl_Tbl_AtomFragment.h	Fl_Tbl_AtomFragment.cxx		\
	Fl_Tbl_Fragment.h		Fl_Tbl_Fragment.cxx			\
	PdfKeyword.h			PdfKeyword.cpp				\
	ProteinDF.h				ProteinDF.cxx				\
	TlOrbitalInfoObject.h	TlOrbitalInfoObject.cpp		\
	TlOrbitalInfo.h			TlOrbitalInfo.cpp			\
	TlOrbitalInfo_Density.h	TlOrbitalInfo_Density.cpp	\
	TlOrbitalInfo_XC.h		TlOrbitalInfo_XC.cpp		\
	TlCombineDensityMatrix.h							\
	TlCombineDensityMatrix.cpp							\
	PdfUserInput.h			PdfUserInput.cpp			\
	PdfUtils.h				PdfUtils.cpp				\
	DfOverlapEngine.h		DfOverlapEngine.cpp			\
	DfHpqEngine.h			DfHpqEngine.cpp				\
	DfEriEngine.h			DfEriEngine.cpp				\
	TlAngularMomentumVector.h							\
	TlAngularMomentumVector.cpp							\
	TlAngularMomentumVectorSet.h						\
	TlAngularMomentumVectorSet.cpp						\
	DfForce.h				DfForce.cpp					\
	DfHpqX.h				DfHpqX.cpp					\
	DfOverlapX.h			DfOverlapX.cpp				\
	DfEriX.h				DfEriX.cpp					\
	TlMatrixCache.h			TlMatrixCache.cpp			\
	DfTaskCtrl.h			DfTaskCtrl.cpp			


libppdf_la_SOURCES = \
	DfInputdata_Parallel.h				DfInputdata_Parallel.cpp			\
	DfIntegrals_Parallel.h				DfIntegrals_Parallel.cpp			\
	DfHpqX_Parallel.h					DfHpqX_Parallel.cpp					\
	DfOverlapX_Parallel.h				DfOverlapX_Parallel.cpp				\
	DfEriX_Parallel.h					DfEriX_Parallel.cpp					\
	DfCD_Parallel.h						DfCD_Parallel.cpp					\
	DfInitialGuess_Parallel.h			DfInitialGuess_Parallel.cpp			\
	DfInitialGuessHarris_Parallel.h		DfInitialGuessHarris_Parallel.cpp	\
	DfInitialGuessHuckel_Parallel.h     DfInitialGuessHuckel_Parallel.cpp	\
	DfPreScf_Parallel.h					DfPreScf_Parallel.cpp				\
	DfScf_Parallel.h					DfScf_Parallel.cpp					\
	DfInvMatrix_Parallel.h				DfInvMatrix_Parallel.cpp			\
	DfGenerateGrid_Parallel.h			DfGenerateGrid_Parallel.cpp			\
	DfDiffDensityMatrix_Parallel.h		DfDiffDensityMatrix_Parallel.cpp	\
	DfDensityFittingX_Parallel.h		DfDensityFittingX_Parallel.cpp		\
	DfDensityFittingX_ScaLAPACK.h		DfDensityFittingX_ScaLAPACK.cpp		\
	DfFockMatrix_Parallel.h				DfFockMatrix_Parallel.cpp			\
	DfGridFreeXC_Parallel.h				DfGridFreeXC_Parallel.cpp			\
	DfXCFunctional_Parallel.h			DfXCFunctional_Parallel.cpp			\
	DfCalcGridX_Parallel.h				DfCalcGridX_Parallel.cpp			\
	DfConvcheck_Parallel.h				DfConvcheck_Parallel.cpp			\
	DfTotalEnergy_Parallel.h			DfTotalEnergy_Parallel.cpp			\
	DfDiagonal_Parallel.h				DfDiagonal_Parallel.cpp				\
	DfXMatrix_Parallel.h				DfXMatrix_Parallel.cpp				\
	DfJMatrix_Parallel.h				DfJMatrix_Parallel.cpp				\
	DfKMatrix_Parallel.h				DfKMatrix_Parallel.cpp				\
	DfTransFmatrix_Parallel.h			DfTransFmatrix_Parallel.cpp			\
	DfTransatob_Parallel.h				DfTransatob_Parallel.cpp			\
	DfDmatrix_Parallel.h				DfDmatrix_Parallel.cpp				\
	DfConverge_Damping_Parallel.h		DfConverge_Damping_Parallel.cpp 	\
	DfConverge_Anderson_Parallel.h		DfConverge_Anderson_Parallel.cpp	\
	DfPopulation_Parallel.h				DfPopulation_Parallel.cpp			\
	DfSummary_Parallel.h				DfSummary_Parallel.cpp				\
	DfCleanup_Parallel.h				DfCleanup_Parallel.cpp				\
	DfForce_Parallel.h					DfForce_Parallel.cpp				\
	DfTaskCtrl_Parallel.h				DfTaskCtrl_Parallel.cpp				\
	ProteinDF_Parallel.h				ProteinDF_Parallel.cpp

if PARALLEL
bin_PROGRAMS = PDF.x PPDF.x
else
bin_PROGRAMS = PDF.x
EXTRA_PROGRAMS = PPDF.x
endif

PDF_x_SOURCES = main.cpp
PDF_x_LDADD = \
	libpdf.la \
	../pdflib/libpdftl.la \
	@LAPACK_LIBS@ @BLAS_LIBS@ @FLIBS@

PPDF_x_SOURCES = main_Parallel.cpp
PPDF_x_LDADD = \
	libpdf.la libppdf.la \
	../pdflib/libppdftl.la \
	../pdflib/libpdftl.la \
	@SCALAPACK_LIBS@ @LAPACK_LIBS@ @BLAS_LIBS@ @MPILIBS@ @FLIBS@

install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	for PROG in $(bin_PROGRAMS); do \
		mv -f $${PROG}$(EXEEXT) $${PROG}-$(VERSION)$(EXEEXT) && \
		$(LN_S) $${PROG}-$(VERSION)$(EXEEXT) $${PROG}$(EXEEXT); \
	done

