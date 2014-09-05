
AC_DEFUN([ACX_SCALAPACK], [
AC_REQUIRE([ACX_BLAS])
AC_REQUIRE([ACX_LAPACK])
AC_REQUIRE([ACX_MPI])
acx_scalapack_ok=no

AC_ARG_WITH(scalapack,
        [AC_HELP_STRING([--with-scalapack=<lib>], [use SCALAPACK library <lib>])])
case $with_scalapack in
        yes | "") ;;
        no) acx_scalapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) SCALAPACK_LIBS="$with_scalapack" ;;
        *) SCALAPACK_LIBS="-l$with_scalapack" ;;
esac

dnl Get fortran linker name of SCALAPACK function to check for.
AC_F77_FUNC(sl_init)

dnl First, check SCALAPACK_LIBS environment variable
if test "x$SCALAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$SCALAPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $MPILIBS $LIBS $FLIBS"
        if test "x$MPICXX" != x; then
                save_CXX="$CXX"; CXX="$MPICXX"
        fi
        AC_MSG_CHECKING([for $sl_init in $SCALAPACK_LIBS])
        AC_TRY_LINK_FUNC($sl_init, [acx_scalapack_ok=yes], [SCALAPACK_LIBS=""])
        AC_MSG_RESULT($acx_scalapack_ok)
        LIBS="$save_LIBS"
        if test "x$MPICXX" != x; then
                CXX="$save_CXX"
        fi
        if test $acx_scalapack_ok = no; then
                SCALAPACK_LIBS=""
        fi
fi


dnl SCALAPACK linked to by default?
if test $acx_scalapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_CHECK_FUNC($sl_init, [acx_scalapack_ok=yes])
        LIBS="$save_LIBS"
fi

dnl Generic SCALAPACK library?
if test $acx_scalapack_ok = no; then
for scalapack in scalapack scalapack_rs6k; do
        if test $acx_scalapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
                AC_CHECK_LIB($scalapack, $sl_init,
                    [acx_scalapack_ok=yes; SCALAPACK_LIBS="-l$scalapack"], [], [$FLIBS])
                LIBS="$save_LIBS"
        fi
done
fi

AC_SUBST(SCALAPACK_LIBS)

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_scalapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_SCALAPACK,1,[Define if you have SCALAPACK library.]),[$1])
        :
else
        acx_scalapack_ok=no
        $2
fi
])dnl ACX_SCALAPACK
