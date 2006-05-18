# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([mincfft],[1.1.2],[Andrew Janke <a.janke@gmail.com>])
AC_CONFIG_SRCDIR([mincfft.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH


# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# Checks for libraries.  See m4/README.
mni_REQUIRE_VOLUMEIO

# check for fftw
AC_CHECK_HEADER(fftw3.h,
   AC_CHECK_LIB([fftw3], fftw_plan_dft_1d, ,
   AC_MSG_ERROR([*** Could not find libfftw3])),
   AC_MSG_ERROR([*** Header file fftw3.h not found])
   )

# Finish up by writing output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT(epm-header)
