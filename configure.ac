dnl Process this file with autoconf to produce a configure script.
AC_INIT(DiagHam,0.01)
AC_CONFIG_SRCDIR([src/config.h])
AM_INIT_AUTOMAKE

AM_CONFIG_HEADER(src/config_ac.h:src/config_ac.h.in)

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AC_PROG_INSTALL

dnl Checks for libraries.
AC_CHECK_LIB(m, cos,MATH_LIBS="-lm",AC_MSG_WARN(math lib not found))
# check if erfc is implemented
AC_CHECK_FUNC(erfc) 

AC_CHECK_LIB(pthread, pthread_create,THREAD_LIBS="-lpthread",AC_MSG_WARN(pthread lib not found : can't use SMP))

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(sys/time.h)

dnl allow for some environment variables to be cached
AC_ARG_VAR(LDFLAGS, "linker flags, e.g. -L<lib dir> if you have libraries in a nonstandard directory <lib dir> (cached"))
AC_ARG_VAR(LIBS, "libraries to pass to the linker, e.g. -l<library> (cached)")
AC_ARG_VAR(CPPFLAGS,"(Objective) C/C++ preprocessor flags, e.g. -I<include dir> if you have headers in a nonstandard directory <include dir> (cached)")
AC_ARG_VAR(CXX,"C++ compiler command (cached)")
AC_ARG_VAR(CXXFLAGS,"C++ compiler flags, help-string (cached)") 

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

dnl parse lapack option

AC_ARG_ENABLE(lapack, [AC_HELP_STRING([--enable-lapack],[enable LAPACK support])], [HAVE_PKG_LAPACK=$enableval],[HAVE_PKG_LAPACK=no])

USER_DEFINED_BLAS=no
AC_ARG_WITH(blas-libs,[AC_HELP_STRING([--with-blas-libs="LIBS"],[BLAS libraries @<:@"-lblas"@:>@])],
[BLAS_LIBS=${withval}
AC_MSG_CHECKING(user-defined BLAS libraries)
AC_MSG_RESULT([${BLAS_LIBS}])
USER_DEFINED_BLAS=yes
]
)

USER_DEFINED_LAPACK=no
AC_ARG_WITH(lapack-libs,[AC_HELP_STRING([--with-lapack-libs="LIBS"],[LAPACK libraries @<:@"-llapack -lg2c"@:>@])],
[LAPACK_LIBS=${withval}
AC_MSG_CHECKING(user-defined LAPACK libraries)
AC_MSG_RESULT([${LAPACK_LIBS}])
USER_DEFINED_LAPACK=yes
]
)

USER_DEFINED_FORTRAN2C=no
AC_ARG_WITH(fortran2c-libs,[AC_HELP_STRING([--with-fortran2c-libs="LIBS"],[FORTRAN2C libraries @<:@"-lgfortran"@:>@])],
[FORTRAN2C_LIBS=${withval}
AC_MSG_CHECKING(user-defined FORTRAN2C libraries)
AC_MSG_RESULT([${FORTRAN2C_LIBS}])
USER_DEFINED_FORTRAN2C=yes
]
)

AC_ARG_WITH(blas-libdir,
[AC_HELP_STRING([--with-blas-libdir=DIR],[BLAS library directory @<:@BLASROOT/lib@:>@])],
[BLAS_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined BLAS library directory)
AC_MSG_RESULT([${BLAS_LIBDIR}])
]
)

AC_ARG_WITH(lapack-libdir,
[AC_HELP_STRING([--with-lapack-libdir=DIR],[LAPACK library directory @<:@LAPACKROOT/lib@:>@])],
[
LAPACK_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined LAPACK library directory)
AC_MSG_RESULT([${LAPACK_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using LAPACK)
AC_MSG_RESULT([${HAVE_PKG_LAPACK}])
 
if test X${HAVE_PKG_LAPACK} = Xyes; then
  AC_DEFINE(HAVE_LAPACK,,[define if we are using LAPACK])
  AC_LANG_PUSH(Fortran 77)
  if test X${USER_DEFINED_BLAS} = Xno; then
	  AC_CHECK_LIB(blas, dgemm, [HAVE_PKG_BLAS=$enableval; BLAS_LIBS="-lblas"],AC_MSG_ERROR(BLAS lib not found))
  fi
  if test X${USER_DEFINED_LAPACK} = Xno; then
     AC_CHECK_LIB(lapack, dsyev, [HAVE_PKG_LAPACK=$enableval; LAPACK_LIBS="-llapack"], AC_MSG_ERROR(LAPACK lib not found), "-lblas")
  fi
  AC_LANG_POP(Fortran 77)
  if test X${USER_DEFINED_FORTRAN2C} = Xno; then
	AC_CHECK_LIB(f2c, main,FORTRAN2C_LIBS="-lf2c",
	AC_CHECK_LIB(g2c, main,FORTRAN2C_LIBS="-lg2c",
	AC_CHECK_LIB(gfortran, main,FORTRAN2C_LIBS="-lgfortran",AC_MSG_ERROR(can't link fortran to c++: no f2c, g2c or gfortran lib found))))
  fi
fi

AM_CONDITIONAL(BUILD_LAPACK, test X${HAVE_PKG_LAPACK} = Xyes)



AC_ARG_ENABLE(lapack-only, [AC_HELP_STRING([--enable-lapack-only],[use LAPACK routines instead of DiagHam ones])], [USE_PKG_LAPACK_ONLY=$enableval],[USE_PKG_LAPACK_ONLY=no])
AC_MSG_CHECKING(whether we are using only LAPACK routines)
AC_MSG_RESULT([${USE_PKG_LAPACK_ONLY}])
 
if test X${USE_PKG_LAPACK_ONLY} = Xyes; then
  AC_DEFINE(HAVE_LAPACK_ONLY,,[define if we are using LAPACK only])
fi


dnl Cheks for gsl library

AC_ARG_ENABLE(gsl, [AC_HELP_STRING([--enable-gsl],[enable gsl support])], [HAVE_PKG_GSL=$enableval],[HAVE_PKG_GSL=no])

AC_MSG_CHECKING(whether we are using gsl library)
AC_MSG_RESULT([${HAVE_PKG_GSL}])


USER_DEFINED_GSL=no
AC_ARG_WITH(gsl-libs,[AC_HELP_STRING([--with-gsl-libs="LIBS"],[GSL libraries @<:@"-lgsl -lgslcblas"@:>@])],
[GSL_LIBS=${withval}
USER_DEFINED_GSL=yes
AC_MSG_CHECKING(user-defined GSL libraries)
AC_MSG_RESULT([${GSL_LIBS}])
])



if test X${USER_DEFINED_GSL} = Xno && test X${HAVE_PKG_GSL} = Xyes; then
  AC_DEFINE(HAVE_GSL,,[define if we are using gsl library])
  AC_CHECK_LIB(gsl,gsl_blas_dgemm,GSL_LIBS="-lgsl -lgslcblas",AC_MSG_ERROR([gsl library not found]),"-lgslcblas")
  GSL_LIBS="-lgsl -lgslcblas"
fi


dnl Cheks for bz2 library

AC_ARG_ENABLE(bz2, [AC_HELP_STRING([--enable-bz2],[enable bz2 support])], [HAVE_PKG_BZ2=$enableval],[HAVE_PKG_BZ2=no])

AC_MSG_CHECKING(whether we are using bz2 library)
AC_MSG_RESULT([${HAVE_PKG_BZ2}])

if test X${HAVE_PKG_BZ2} = Xyes; then
  AC_DEFINE(HAVE_BZ2,,[define if we are using bzip2 library])
  AC_CHECK_LIB(bz2,BZ2_bzReadOpen,BZ2_LIBS=-lbz2,AC_MSG_ERROR([bzip2 library not found]))
fi

dnl Cheks for gmp library

AC_ARG_ENABLE(gmp, [AC_HELP_STRING([--enable-gmp],[enable gmp support])], [HAVE_PKG_GMP=$enableval],[HAVE_PKG_GMP=no])

AC_MSG_CHECKING(whether we are using gmp library)
AC_MSG_RESULT([${HAVE_PKG_GMP}])


USER_DEFINED_GMP=no
AC_ARG_WITH(gmp-libs,[AC_HELP_STRING([--with-gmp-libs="LIBS"],[GMP libraries @<:@"-lgmp -lgmpxx"@:>@])],
[GMP_LIBS=${withval}
USER_DEFINED_GMP=yes
AC_MSG_CHECKING(user-defined GMP libraries)
AC_MSG_RESULT([${GMP_LIBS}])
])



if test X${USER_DEFINED_GMP} = Xno && test X${HAVE_PKG_GMP} = Xyes; then
  AC_DEFINE(HAVE_GMP,,[define if we are using gmp library])
  AC_CHECK_LIB(gmp,__gmpz_init,GMP_LIBS=-lgmp,AC_MSG_ERROR([gmp library not found]))
  GMP_LIBS="-lgmp -lgmpxx"
fi


dnl Cheks for mpack arbitrary precision LAPACK library

AC_ARG_ENABLE(mpack, [AC_HELP_STRING([--enable-mpack],[enable mpack support])], [HAVE_PKG_MPACK=$enableval],[HAVE_PKG_MPACK=no])

AC_MSG_CHECKING(whether we are using mpack library)
AC_MSG_RESULT([${HAVE_PKG_MPACK}])


USER_DEFINED_MPACK=no
AC_ARG_WITH(mpack-libs,[AC_HELP_STRING([--with-mpack-libs="LIBS"],[MPACK libraries @<:@"-lmlapack_gmp -lmblas_gmp"@:>@])],
[MPACK_LIBS=${withval}
USER_DEFINED_MPACK=yes
AC_MSG_CHECKING(user-defined MPACK libraries)
AC_MSG_RESULT([${MPACK_LIBS}])
])

if test X${USER_DEFINED_MPACK} = Xno && test X${HAVE_PKG_MPACK} = Xyes; then
  AC_DEFINE(HAVE_MPACK,,[define if we are using mpack library])
#  AC_CHECK_LIB(mlapack_gmp,__Chpev,MPACK_LIBS=-lmlapack_gmp,AC_MSG_ERROR([mpack library not found])) # edit this line for a proper check
  MPACK_LIBS="-lmlapack_gmp -lmblas_gmp"
fi

if test X${HAVE_PKG_MPACK} = Xyes; then
  AC_MSG_CHECKING(augmented LAPACK linker commands with MPACK)
  LAPACK_LIBS="${LAPACK_LIBS} ${MPACK_LIBS}"
  AC_MSG_RESULT([${LAPACK_LIBS}])
fi

dnl Cheks for fftw3 library

AC_ARG_ENABLE(fftw, [AC_HELP_STRING([--enable-fftw],[enable fftw support])], [HAVE_PKG_FFTW=$enableval],[HAVE_PKG_FFTW=no])

AC_MSG_CHECKING(whether we are using fftw library)
AC_MSG_RESULT([${HAVE_PKG_FFTW}])

if test X${HAVE_PKG_FFTW} = Xyes; then
  AC_DEFINE(HAVE_FFTW,,[define if we are using fftw library])
  AC_CHECK_LIB(fftw3,fftw_plan_dft_r2c_1d,FFTW_LIBS=-lfftw3,AC_MSG_ERROR([fftw library not found]))
  FFTW_LIBS="-lfftw3"
fi


dnl parse mpi options

AC_ARG_ENABLE(mpi, [AC_HELP_STRING([--enable-mpi],[enable MPI support])], [HAVE_PKG_MPI=$enableval],[HAVE_PKG_MPI=no])
AC_ARG_WITH(mpi-cxx, [AC_HELP_STRING([--with-mpi-cxx], [use MPI C++ compiler (enables MPI) @<:@mpiCC@:>@])],
[HAVE_PKG_MPI=yes 
if test X${withval} = Xyes; then
  MPI_CXX=mpiCC
else
  MPI_CXX=${withval}
fi
AC_MSG_CHECKING(user-defined MPI C++ compiler)
AC_MSG_RESULT([${MPI_CXX}])
]
)

AC_ARG_WITH(mpi-libs,[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[MPI_LIBS=${withval}
AC_MSG_CHECKING(user-defined MPI libraries)
AC_MSG_RESULT([${MPI_LIBS}])
]
)


AC_ARG_WITH(mpi-incdir,[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@])],
[
MPI_INC=${withval}
AC_MSG_CHECKING(user-defined MPI includes)
AC_MSG_RESULT([${MPI_INC}])
]
)


AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@])],
[
MPI_LIBDIR=${withval}
AC_MSG_CHECKING(user-defined MPI library directory)
AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])
 
if test X${HAVE_PKG_MPI} = Xyes; then
  AC_DEFINE(HAVE_MPI,,[define if we are using MPI])
fi

AM_CONDITIONAL(BUILD_MPI, test X${HAVE_PKG_MPI} = Xyes)

dnl check for MPI c++ compiler

if test X${HAVE_PKG_MPI} = Xyes; then

  if test -z "${MPI_DIR}" &&  test -z "${MPI_CXX}"; then
    MPI_CXX=mpiCC
  fi
fi

if test -n "${MPI_CXX}"; then

  if test -f ${MPI_CXX}; then   
    MPI_CXX_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  fi

  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi

fi

dnl parse scalapack option

#AC_ARG_ENABLE(scalapack, [AC_HELP_STRING([--enable-scalapack],[enable SCALAPACK support])], [HAVE_PKG_SCALAPACK=$enableval],[HAVE_PKG_SCALAPACK=no])
#
#USER_DEFINED_SCALAPACK=no
#AC_ARG_WITH(scalapack-libs,[AC_HELP_STRING([--with-scalapack-libs="LIBS"],[SCALAPACK libraries @<:@"-lscalapack "@:>@])],
#[SCALAPACK_LIBS=${withval}
#AC_MSG_CHECKING(user-defined SCALAPACK libraries)
#AC_MSG_RESULT([${SCALAPACK_LIBS}])
#USER_DEFINED_SCALAPACK=yes
#]
#)
#
#AC_MSG_CHECKING(whether we are using SCALAPACK)
#AC_MSG_RESULT([${HAVE_PKG_SCALAPACK}])
#
#if test X${HAVE_PKG_SCALAPACK} = Xyes; then
# AC_DEFINE(HAVE_SCALAPACK,,[define if we are using SCALAPACK])
# AC_CHECK_LIB(scalapack, pdsytrd, [HAVE_PKG_SCALAPACK=$enableval; SCALAPACK_LIBS="-lscalapack"], AC_MSG_ERROR(SCALAPACK lib not found))
#fi


AC_ARG_ENABLE(scalapack, [AC_HELP_STRING([--enable-scalapack],[enable SCALAPACK support])], [HAVE_PKG_SCALAPACK=$enableval],[HAVE_PKG_SCALAPACK=no])

USER_DEFINED_SCALAPACK=no
AC_ARG_WITH(scalapack-libs,[AC_HELP_STRING([--with-scalapack-libs="LIBS"],[SCALAPACK libraries @<:@"-lscalapack -lg2c"@:>@])],
[SCALAPACK_LIBS=${withval}
AC_MSG_CHECKING(user-defined SCALAPACK libraries)
AC_MSG_RESULT([${SCALAPACK_LIBS}])
USER_DEFINED_SCALAPACK=yes
]
)

AC_ARG_WITH(scalapack-libdir,
[AC_HELP_STRING([--with-scalapack-libdir=DIR],[SCALAPACK library directory @<:@SCALAPACKROOT/lib@:>@])],
[
SCALAPACK_LIBDIR=-L${withval}
AC_MSG_CHECKING(user-defined SCALAPACK library directory)
AC_MSG_RESULT([${SCALAPACK_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using SCALAPACK)
AC_MSG_RESULT([${HAVE_PKG_SCALAPACK}])

if test X${HAVE_PKG_SCALAPACK} = Xyes; then
  if test X${HAVE_PKG_MPI} = Xno; then
    AC_MSG_ERROR([SCALAPACK requires MPI.])
  fi
  if test X${HAVE_PKG_LAPACK} = Xno; then
    AC_MSG_ERROR([SCALAPACK requires LAPACK.])
  fi
fi

if test X${HAVE_PKG_SCALAPACK} = Xyes; then
 AC_DEFINE(HAVE_SCALAPACK,,[define if we are using SCALAPACK])
 AC_LANG_PUSH([Fortran])
 if test X${USER_DEFINED_SCALAPACK} = Xno; then
    AC_CHECK_LIB(scalapack, pdsytrd, [HAVE_PKG_SCALAPACK=$enableval; SCALAPACK_LIBS="-lscalapack"], AC_MSG_ERROR(SCALAPACK lib not found), ${SCALAPACK_LIBDIR} -lblas -llapack)
 fi
 AC_LANG_POP([Fortran])
 if test X${USER_DEFINED_LAPACK} = Xno; then
	AC_CHECK_LIB(f2c, main,FORTRAN2C_LIBS="-lf2c",
	AC_CHECK_LIB(g2c, main,FORTRAN2C_LIBS="-lg2c",
	AC_CHECK_LIB(gfortran, main,FORTRAN2C_LIBS="-lgfortran",AC_MSG_ERROR(can't link fortran to c++: no f2c, g2c or gfortran lib found))))
  fi
fi

AM_CONDITIONAL(BUILD_SCALAPACK, test X${HAVE_PKG_SCALAPACK} = Xyes)

dnl check for intel mkl 

AC_ARG_ENABLE(intelmkl, [AC_HELP_STRING([--enable-intelmkl],[using intel-mkl librairies for LAPACK support])], [HAVE_PKG_INTELMKL=$enableval],[HAVE_PKG_INTELMKL=no])

AC_MSG_CHECKING(whether we are using INTEL-MKL)
AC_MSG_RESULT([${HAVE_PKG_INTELMKL}])

 AC_ARG_WITH(intelmkl-libdir,
  [AC_HELP_STRING([--with-intelmkl-libdir=DIR],[INTEL MKL library directory @<:@INTELMKLROOT/lib@:>@])],
  [
  INTELMKL_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined INTELMKL library directory)
  AC_MSG_RESULT([${INTELMKL_LIBDIR}])
  ])
 
if test X${HAVE_PKG_INTELMKL} = Xyes; then
  AC_CHECK_LIB([mkl_intel_ilp64],fftw_execute,INTELMKL_LIBS="-L${INTELMKL_LIBDIR} -lmkl_intel_thread -lmkl_core -lmkl_intel_ilp64 -lpthread -lm -ldl ",
  AC_MSG_ERROR([Cannot find Intel Math Kernel Library >= 6.0]),[-L${INTELMKL_LIBDIR} -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -openmp -liomp5  --lpthread])
  AC_DEFINE(HAVE_LAPACK,,[define if we are using LAPACK])
fi



dnl check byte ordering


ARG_GLOBAL_COMMAND_LOG=none
AC_ARG_WITH(cmd-log,[AC_HELP_STRING([--with-cmd-log="LOGFILE"],[global logfile where DiagHam commands are stored by default])],
[ARG_GLOBAL_COMMAND_LOG=${withval}
AC_MSG_CHECKING(global command log)
AC_MSG_RESULT([${ARG_GLOBAL_COMMAND_LOG}])
]
)

if test X${ARG_GLOBAL_COMMAND_LOG} != Xnone; then
   AC_DEFINE(HAVE_GLOBAL_COMMAND_LOG, [],[using global logfile ])
   AC_DEFINE_UNQUOTED(GLOBAL_COMMAND_LOG, ["${ARG_GLOBAL_COMMAND_LOG}"],[location of global command log])
fi


AC_TRY_RUN([
int main()
{
	short s = 1;
	short* ptr = &s;
	unsigned char c = *((char*)ptr);
	return c;	
}

]
, AC_DEFINE(HAVE_BIGENDIAN, [], [use big endian byte ordering])
, AC_DEFINE(HAVE_LITTLEENDIAN, [], [use little endian byte ordering])
)


dnl Checks for library functions.

AC_MSG_CHECKING(environment for CXXFLAGS)

if test -n "${CXXFLAGS}"; then
  AC_MSG_RESULT([${CXXFLAGS}])
else
  CXXFLAGS=" -O2 -Wall "
  AC_MSG_RESULT([using defaults])
fi


AM_CPPFLAGS="-I\$(top_srcdir)/src -I\$(top_srcdir)/Base/src"

dnl check for debugging/profiling mode 

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[enable debug mode])], [HAVE_DEBUG=$enableval
CXXFLAGS="${CXXFLAGS} -g -fPIC"
],[HAVE_DEBUG=no])

AC_ARG_ENABLE(profile, [AC_HELP_STRING([--enable-profile],[enable profile mode])], [HAVE_PROFILE=$enableval
CXXFLAGS="${CXXFLAGS} -g -pg"
],[HAVE_PROFILE=no])

if test X${HAVE_PKG_INTELMKL} = Xyes; then
CXXFLAGS="${CXXFLAGS} -DMKL_ILP64 -qopenmp -I${MKLROOT}/include"
fi

AC_MSG_CHECKING(whether we are using debug mode)
AC_MSG_RESULT([${HAVE_DEBUG}])

AC_MSG_CHECKING(whether we are using profile mode)
AC_MSG_RESULT([${HAVE_PROFILE}])

dnl force use of -m64 option

AC_ARG_ENABLE(m64, [AC_HELP_STRING([--enable-m64],[force -m64 mode])], [HAVE_FORCEM64=$enableval
CXXFLAGS="${CXXFLAGS} -m64"
],[HAVE_FORCEM64=no])

AC_MSG_CHECKING(whether we are forcing -m64 mode)
AC_MSG_RESULT([${HAVE_FORCEM64}])


dnl check which packages to compile

AC_ARG_ENABLE(quantumdots, [AC_HELP_STRING([--enable-quantumdots],[enable the quantum dot package])], [HAVE_QUANTUMDOTS=yes
AM_CPPFLAGS="${AM_CPPFLAGS} -I\$(top_srcdir)/QuantumDots/src"
], [HAVE_QUANTUMDOTS=no])

AC_ARG_ENABLE(fqhe, [AC_HELP_STRING([--enable-fqhe],[enable the fractional quantum Hall effect package])], [HAVE_FQHE=yes
AM_CPPFLAGS="${AM_CPPFLAGS} -I\$(top_srcdir)/FQHE/src"
], [HAVE_FQHE=no])

AC_ARG_ENABLE(fti, [AC_HELP_STRING([--enable-fti],[enable the fractional topological insulator package])], [HAVE_FTI=yes
AM_CPPFLAGS="${AM_CPPFLAGS} -I\$(top_srcdir)/FTI/src"
AC_DEFINE(HAVE_FTI, [], [have fractional topological order package])
], [HAVE_FTI=no])

AC_ARG_ENABLE(spin, [AC_HELP_STRING([--enable-spin],[enable the spin system package])], [HAVE_SPIN=yes
AM_CPPFLAGS="${AM_CPPFLAGS} -I\$(top_srcdir)/Spin/src"
], [HAVE_SPIN=no])

AC_ARG_ENABLE(anyons, [AC_HELP_STRING([--enable-anyons],[enable the anyonic model package])], [HAVE_ANYONS=yes
AM_CPPFLAGS="${AM_CPPFLAGS} -I\$(top_srcdir)/Anyons/src"
], [HAVE_ANYONS=no])

AC_MSG_CHECKING(whether we are compiling the quantum dot package)
AC_MSG_RESULT([${HAVE_QUANTUMDOTS}])

AC_MSG_CHECKING(whether we are compiling the spin system package)
AC_MSG_RESULT([${HAVE_SPIN}])

AC_MSG_CHECKING(whether we are compiling the fractional quantum Hall effect package)
AC_MSG_RESULT([${HAVE_FQHE}])

AC_MSG_CHECKING(whether we are compiling the fractional topological insulator package)
AC_MSG_RESULT([${HAVE_FTI}])

if test "$HAVE_QUANTUMDOTS" = yes; then
       MAYBE_QUANTUMDOTS="QuantumDots"
     else
       MAYBE_QUANTUMDOTS=""
fi

if test "$HAVE_SPIN" = yes; then
       MAYBE_SPIN="Spin"
     else
       MAYBE_SPIN=""
fi

if test "$HAVE_FQHE" = yes; then
       MAYBE_FQHE="FQHE"
     else
       MAYBE_FQHE=""
fi

if test "$HAVE_FTI" = yes; then
       MAYBE_FTI="FTI"
     else
       MAYBE_FTI=""
fi

if test "$HAVE_ANYONS" = yes; then
       MAYBE_ANYONS="Anyons"
     else
       MAYBE_ANYONS=""
fi


dnl add further includes for GSL library if needed

if test X${HAVE_PKG_GSL} = Xyes; then
  AM_CPPFLAGS="${GSL_AM_CPPFLAGS} ${AM_CPPFLAGS}"
fi


#dnl add further includes for ARPACK library if needed
#
#if test X${HAVE_PKG_ARPACK} = Xyes; then
#  AM_CPPFLAGS="${AM_CPPFLAGS} -I${ARPACKPP_INCLUDE}"
#  if test X${USER_DEFINED_ARPACK_DIR} = Xyes; then
#    ARPACK_LIBS="-L${ARPACK_LIBDIR} ${ARPACK_LIBS}"
#  fi
#fi	


dnl end of configure

AC_SUBST(CXXFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(MATH_LIBS)
AC_SUBST(THREAD_LIBS)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(BLAS_LIBDIR)
AC_SUBST(LAPACK_LIBDIR)
AC_SUBST(SCALAPACK_LIBS)
AC_SUBST(SCALAPACK_LIBDIR)
#AC_SUBST(ARPACK_LIBS)
AC_SUBST(FORTRAN2C_LIBS)
AC_SUBST(BZ2_LIBS)
AC_SUBST(GSL_LIBS)
AC_SUBST(GMP_LIBS)
AC_SUBST(MPACK_LIBS)
AC_SUBST(FFTW_LIBS)
AC_SUBST(INTELMKL_LIBS)
AC_SUBST(MAYBE_QUANTUMDOTS)
AC_SUBST(MAYBE_SPIN)
AC_SUBST(MAYBE_FQHE)
AC_SUBST(MAYBE_FTI)
#AC_SUBST(MAYBE_ANYONS)

AC_CONFIG_FILES(Makefile src/Makefile src/Programs/Makefile src/Architecture/Makefile src/Architecture/ArchitectureOperation/Makefile src/FunctionBasis/Makefile src/Tensor/Makefile src/QuantumNumber/Makefile src/Interaction/Makefile src/Interaction/InternalInteraction/Makefile src/HilbertSpace/Makefile src/HilbertSpace/DMRGHilbertSpace/Makefile src/HilbertSpace/ManyBodyHilbertSpace/Makefile src/Output/Makefile src/GeneralTools/Makefile src/Vector/Makefile src/Polynomial/Makefile src/Matrix/Makefile src/MCObservables/Makefile src/LanczosAlgorithm/Makefile src/Hamiltonian/Makefile src/Hamiltonian/DMRGHamiltonian/Makefile src/TensorProduct/Makefile src/DMRGAlgorithm/Makefile src/Operator/Makefile src/MathTools/Makefile src/MathTools/NumericalAnalysis/Makefile src/MathTools/RandomNumber/Makefile src/Options/Makefile src/Tools/Makefile src/MainTask/Makefile src/Architecture/ClusterArchitecture/Makefile src/MPSObjects/Makefile Base/Makefile Base/src/Makefile Base/src/BitmapTools/Makefile Base/src/BitmapTools/Color/Makefile Base/src/BitmapTools/BitmapPicture/Makefile QuantumDots/Makefile QuantumDots/src/Makefile QuantumDots/src/Tools/Makefile QuantumDots/src/HilbertSpace/Makefile QuantumDots/src/Hamiltonian/Makefile QuantumDots/src/Tools/Spectra/Makefile QuantumDots/src/Tools/Potential/Makefile QuantumDots/src/Tools/Analysis/Makefile QuantumDots/src/Programs/Makefile FQHE/Makefile FQHE/src/Makefile FQHE/src/Architecture/Makefile FQHE/src/Tools/Makefile FQHE/src/Programs/Makefile FQHE/src/Programs/FQHEOnSphere/Makefile FQHE/src/Programs/FQHEOnDisk/Makefile FQHE/src/Programs/FQHEOnTorus/Makefile FQHE/src/Programs/FQHEOnCylinder/Makefile FQHE/src/Programs/FQHEOnLattice/Makefile FQHE/src/Architecture/ArchitectureOperation/Makefile FQHE/src/FunctionBasis/Makefile FQHE/src/QuantumNumber/Makefile FQHE/src/HilbertSpace/Makefile FQHE/src/Hamiltonian/Makefile FQHE/src/MainTask/Makefile FQHE/src/Tools/FQHEWaveFunction/Makefile FQHE/src/Tools/FQHEMonteCarlo/Makefile FQHE/src/Tools/FQHEMPS/Makefile FQHE/src/Tools/FQHESpectrum/Makefile FQHE/src/Tools/FQHEFiles/Makefile FQHE/src/Operator/Makefile Spin/Makefile Spin/src/Makefile Spin/src/MPSObjects/Makefile Spin/src/Architecture/Makefile Spin/src/Programs/Makefile Spin/src/Architecture/ArchitectureOperation/Makefile Spin/src/HilbertSpace/Makefile Spin/src/Hamiltonian/Makefile Spin/src/Operator/Makefile Spin/src/Tools/Makefile Spin/src/Tools/SpinFiles/Makefile FTI/Makefile FTI/src/Makefile FTI/src/Architecture/Makefile FTI/src/Architecture/ArchitectureOperation/Makefile FTI/src/HilbertSpace/Makefile FTI/src/Hamiltonian/Makefile FTI/src/Tools/Makefile FTI/src/Operator/Makefile FTI/src/Tools/FTITightBinding/Makefile FTI/src/Tools/FTIFiles/Makefile FTI/src/Programs/Makefile FTI/src/Programs/FCI/Makefile FTI/src/Programs/FTI/Makefile FTI/src/Programs/HubbardModels/Makefile)
# Anyons/Makefile Anyons/src/Makefile Anyons/src/HilbertSpace/Makefile Anyons/src/Programs/Makefile)

AC_OUTPUT
