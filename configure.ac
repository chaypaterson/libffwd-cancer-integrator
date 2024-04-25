AC_INIT([cancer-integrator],[0.1],[chay.paterson@manchester.ac.uk])
#AM_INIT_AUTOMAKE([foreign subdir-objects])

#TODO dependencies: GSL, Eigen
# TODO: will need to ask for /include/ dirs for eigen
# should test for eigen version

# Check for dependencies:
# Check for GSL functions:
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

AM_INIT_AUTOMAKE
AC_PROG_CXX
AC_CONFIG_FILES([Makefile])
AC_ENABLE_SHARED
LT_INIT
AC_OUTPUT
