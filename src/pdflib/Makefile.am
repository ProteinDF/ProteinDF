AM_CPPFLAGS= -I$(top_srcdir)/include 

libpdftl_la_SOURCES = \
	TlFname.cxx 					\
	TlSystem.cpp					\
	TlFile.cpp						\
	TlStringTokenizer.cpp			\
	TlUtils.cpp						\
	TlSerializeData.cpp				\
	TlPseudoYaml.cpp				\
	TlMsgPack.cpp					\
	TlTime.cpp						\
	TlParameter.cpp					\
	TlMath.cxx						\
	TlPosition.cpp	TlAtom.cpp 		\
	TlPrdctbl.cxx	TlResidue.cxx	\
	TlMemManager.cpp				\
	\
	TlMatrixObject.cpp				\
	TlMatrix.cpp 					\
	TlSparseMatrix.cpp				\
	TlSparseSymmetricMatrix.cpp		\
	TlSparseVectorMatrix.cpp		\
	TlSymmetricMatrix.cpp			\
	TlPartialMatrix.cpp				\
	TlPartialSymmetricMatrix.cpp	\
	TlFileMatrix.cpp				\
	TlFileSymmetricMatrix.cpp		\
	TlListMatrix.cpp				\
	TlRowVectorMatrix.cpp			\
	TlColVectorMatrix.cpp			\
	TlRowVectorMatrix2.cpp 			\
	TlColVectorMatrix2.cpp 			\
	TlVector.cpp					\
	TlFmt.cpp						\
	TlLogging.cpp					\
	TlGetopt.cpp

libppdftl_la_SOURCES = \
	TlCommunicate.cpp				\
	TlDistributeMatrix.cpp			\
	TlDistributeSymmetricMatrix.cpp	\
	TlDistributeVector.cpp		


if PARALLEL
lib_LTLIBRARIES = libpdftl.la libppdftl.la
else
lib_LTLIBRARIES = libpdftl.la
EXTRA_LTLIBRARIES = libppdftl.la
endif


##	TlSparseHashMatrix.cpp		
##	TlSparseHashSymmetricMatrix.cpp	
