AC_INIT([cancer-integrator],[0.8],[chay.paterson@manchester.ac.uk])
#AM_INIT_AUTOMAKE([foreign subdir-objects])

## Set the language to C++:
AC_LANG_PUSH([C++])
: ${CXXFLAGS="-O3 --std=c++17"}

## Check for dependencies:
## Check for GSL functions:
AC_CHECK_LIB([gsl], [gsl_rng_alloc], [], [
    AC_MSG_ERROR([gsl_rng_alloc not found: is GSL installed?])
])
AC_CHECK_LIB([gsl], [gsl_rng_env_setup], [], [
    AC_MSG_ERROR([gsl_rng_env_setup not found: is GSL installed?])
])
AC_CHECK_LIB([gsl], [gsl_ran_exponential], [], [
    AC_MSG_ERROR([gsl_ran_exponential not found: is GSL installed?])
])
AC_CHECK_LIB([gsl], [gsl_ran_flat], [], [
    AC_MSG_ERROR([gsl_ran_flat not found: is GSL installed?])
])
AC_CHECK_LIB([gsl], [gsl_rng_free], [], [
    AC_MSG_ERROR([gsl_rng_free not found: is GSL installed?])
])
AC_CHECK_LIB([gsl], [gsl_rng_set], [], [
    AC_MSG_ERROR([gsl_rng_set not found: is GSL installed?])
])
## Check for GSL headers:
AC_CHECK_HEADERS([gsl/gsl_rng.h gsl/gsl_randist.h], [], [
    AC_MSG_ERROR([Required GSL RNG headers not found: is GSL installed?])
])

## Check for Eigen headers:
AC_ARG_WITH([eigen], AS_HELP_STRING([--with-eigen=...], [path to Eigen 3]),
CPPFLAGS="$CPPFLAGS -I$with_eigen",
[])

## (NB Eigen is a header-only library)
AC_CHECK_HEADERS([Eigen/Dense Eigen/Eigenvalues], [], [
    AC_MSG_ERROR([Required Eigen headers not found: is Eigen 3 installed? Check config.log for details.])
])
AC_LANG_POP([C++])

AM_INIT_AUTOMAKE
AC_PROG_CXX
AC_CONFIG_FILES([Makefile])
AC_ENABLE_SHARED
LT_INIT
AC_OUTPUT
