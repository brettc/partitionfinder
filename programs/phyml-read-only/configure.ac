#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

dnl AC_PREREQ(2.61)
dnl AC_INIT([PhyML(beta)],m4_esyscmd([sh -c "head -n 1 version | sed -e 's|^|v3.0-|g' | sed -e 's|$||g' | tr -d '\n'"]), [guindon@stat.auckland.ac.nz])
AC_INIT([PhyML],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[s.guindon@auckland.ac.nz])
dnl AC_INIT([PhyML],[beta],[guindon@stat.auckland.ac.nz])
dnl AC_CONFIG_SRCDIR([simu.c])
AC_CONFIG_SRCDIR([src/simu.c])
AC_CONFIG_HEADERS([config.h])

AC_DEFINE([UNIX],[1],[Unix tag on])
AC_DEFINE([DEBUG],[1],[Debug tag on])

AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_SETVBUF_REVERSED
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor pow rint sqrt strchr strstr])

# Thanks to Shane Sturrock from Biomatters
#case "${host}" in
#*darwin*)
#AM_INIT_AUTOMAKE([no-dependencies])
ARCH_flag="-arch i386 -mmacosx-version-min=10.4"
#;;
#esac


dnl CFLAGS="-O3 -fomit-frame-pointer -funroll-loops -Wall ${ARCH_flag}"
dnl CFLAGS="-O3 -Wfloat-equal -fomit-frame-pointer -funroll-loops"
dnl CFLAGS="-O3 -Wfloat-equal"
dnl CFLAGS="-g"

LDFLAGS="${ARGC_flag}"

dnl LT_INIT
dnl AC_PROG_LIBTOOL
dnl Add option for compiling with debugging informations

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [Remove optimization options and add debug informations.])])
AS_IF([test "x$enable_debug" = "xyes"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O0 -g"],
      [CFLAGS="-ansi -pedantic -Wall -std=c99 -O3 -fomit-frame-pointer -funroll-loops ${ARCH_flag}"]) 



AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi],
                              [Compile with mpicc instead of gcc.])])
AS_IF([test "x$enable_mpi" = "xyes"],[CC="mpicc"])
AS_IF([test "x$enable_mpi" = "xyes"],AC_DEFINE([MPI],[1],[MPI tag on]))
AM_CONDITIONAL([WANT_MPI], [test "x$enable_mpi" = "xyes"]) 

AC_ARG_ENABLE([phytime],[AS_HELP_STRING([--enable-phytime],[Compile PhyTime])],[phytime=yes],[phytime=no])
AM_CONDITIONAL([WANT_PHYTIME], [test "$phytime" = yes]) 
if test "$phytime" = yes; then
   AC_DEFINE([PHYTIME],[1],[PHYTIME tag on])
fi

AC_ARG_ENABLE([tiporder],[AS_HELP_STRING([--enable-tiporder],[Compile tiporder])],[tiporder=yes],[tiporder=no])
AM_CONDITIONAL([WANT_TIPORDER], [test "$tiporder" = yes]) 
if test "$tiporder" = yes; then
   AC_DEFINE([TIPORDER],[1],[TIPORDER tag on])
fi

AC_ARG_ENABLE([part],[AS_HELP_STRING([--enable-part],[Compile Part])],[part=yes],[part=no])
AM_CONDITIONAL([WANT_PART], [test "$part" = yes]) 
if test "$part" = yes; then
   AC_DEFINE([PART],[1],[PART tag on])
fi

AC_ARG_ENABLE([rwrap],[AS_HELP_STRING([--enable-rwrap],[Compile Rwrap])],[rwrap=yes],[rwrap=no])
AM_CONDITIONAL([WANT_RWRAP], [test "$rwrap" = yes]) 
if test "$rwrap" = yes; then
   AC_DEFINE([RWRAP],[1],[RWRAP tag on])
fi

AC_ARG_ENABLE([phycont],[AS_HELP_STRING([--enable-phycont],[Compile PhyCont])],[phycont=yes],[phycont=no])
AM_CONDITIONAL([WANT_PHYCONT], [test "$phycont" = yes]) 
if test "$phycont" = yes; then
   AC_DEFINE([PHYCONT],[1],[PHYCONT tag on])
fi

AC_ARG_ENABLE([m4],[AS_HELP_STRING([--enable-m4],[Compile M4])],[m4=yes],[m4=no])
AM_CONDITIONAL([WANT_M4], [test "$m4" = yes]) 
if test "$m4" = yes; then
   AC_DEFINE([M4],[1],[m4 tag on])
fi

AC_ARG_ENABLE([rf],[AS_HELP_STRING([--enable-rf],[Compile RF])],[rf=yes],[rf=no])
AM_CONDITIONAL([WANT_RF], [test "$rf" = yes]) 
if test "$rf" = yes; then
   AC_DEFINE([RF],[1],[RF tag on])
fi

if test "$phytime" = no; then
if test "$tiporder" = no; then
if test "$part" = no; then
if test "$rwrap" = no; then
if test "$phycont" = no; then
if test "$m4" = no; then
   AC_DEFINE([PHYML],[1],[PHYML tag on])
fi
fi
fi
fi
fi
fi

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
