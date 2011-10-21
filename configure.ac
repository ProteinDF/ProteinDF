#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ProteinDF], [2011.1git]m4_esyscmd([git log -1 --format="%h" | tr -d "\n"]))
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/pdf/main.cpp])
AC_CONFIG_MACRO_DIR([m4])

# ====================================================================
# config.h
# ====================================================================

# if you want to use SSL2, uncomment below and do bootstrap.sh.
#AC_DEFINE_UNQUOTED(main, MAIN__, [for SSL2])

PACKAGE=$PACKAGE_NAME

# Checks for programs.
AC_LANG([C++])

AC_PROG_F77
AC_FC_LIBRARY_LDFLAGS
AC_FC_MAIN
AC_FC_WRAPPERS

AC_PROG_CXX

AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL
AC_PROG_INSTALL

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

AC_CXX_HAVE_STL
AC_CXX_HAVE_VALARRAY

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_CHECK_FUNCS([floor getcwd pow sqrt])

AC_C_COMPILE_VALUE(sizeof(int), sizeof_int)
AC_C_COMPILE_VALUE(sizeof(long), sizeof_long)
AX_CREATE_STDINT_H(include/pdflib-int.h)

# ====================================================================
# for sleep
# ====================================================================
AX_SHORT_SLEEP([WITH_SLEEP=yes],[WITH_SLEEP=no])

# ====================================================================
# for LAPACK
# ====================================================================
ACX_BLAS
ACX_LAPACK

# ====================================================================
# for MPI
# ====================================================================
ACX_MPI
if test "x$MPICXX" != x; then
   CXX="$MPICXX"
fi

# ====================================================================
# for ScaLAPACK
# ====================================================================
ACX_SCALAPACK

# ====================================================================
# for parallel
# ====================================================================
AC_ARG_ENABLE(parallel, [  --enable-parallel  build PPDF.x],
[case "${enableval}" in
      yes) parallel=true ;;
      no)  parallel=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-parallel) ;;
esac], [parallel=false])
AM_CONDITIONAL(PARALLEL, test x$parallel = xtrue)

# ====================================================================
# for DEBUG
# ====================================================================
AC_ARG_ENABLE(debug, [  --enable-debug  trun on debugging [default no]],,enable_debug=no)
AC_MSG_CHECKING(whether to enable debuging)
if test x$enable_debug = xyes; then
   AC_MSG_RESULT(yes)
   CXXFLAGS="-g ${CXXFLAGS}"
else
   AC_MSG_RESULT(no)
   AC_DEFINE(NDEBUG, [1], [turn off assert])
#   CXXFLAGS=${CXXFLAGS/-g/}
fi
AC_SUBST(CXXFLAGS)

# ====================================================================
# for PROFILE
# ====================================================================
AC_ARG_ENABLE(profile, [  --enable-profile  trun on profiling [default no]],,enable_profile=no)
AC_MSG_CHECKING(whether to enable profiling)
if test x$enable_profile = xyes; then
   AC_MSG_RESULT(yes)
   CXXFLAGS="${CXXFLAGS} -pg"
else
   AC_MSG_RESULT(no)
#   CXXFLAGS=${CXXFLAGS/-pg/}
fi
AC_SUBST(CXXFLAGS)

# ====================================================================
# for CppUnit
# ====================================================================
AC_ARG_ENABLE(cppunit,
[  --enable-cppunit  turn on cppunit test],
[case "${enableval}" in
      yes) cppunit=true ;;
      no)  cppunit=false ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-cppunit) ;;
esac], [cppunit=false])
AM_CONDITIONAL(CPPUNIT, test x$cppunit = xtrue)


AM_CONFIG_HEADER([include/config.h])
AC_CONFIG_FILES([
	Makefile
	include/Makefile
	src/Makefile
	src/pdflib/Makefile
	src/pdf/Makefile
	src/tools/Makefile
	src/unit_test/Makefile
	src/parallel_test/Makefile
	src/performace_test/Makefile
	scripts/Makefile
	data/Makefile
	])
AC_OUTPUT

