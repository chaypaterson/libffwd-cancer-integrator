AC_INIT([cancer-integrator],[0.1],[chay.paterson@manchester.ac.uk])
AM_INIT_AUTOMAKE([foreign subdir-objects])

#TODO dependencies: GSL, Eigen
# TODO: will need to ask for /include/ dirs for eigen
# should test for eigen version

LT_INIT
AC_CONFIG_MACRO_DIRS([m4])
ACLOCAL_AMFLAGS = -I m4

AC_PROG_CXX
AC_CONFIG_FILES([
 Makefile
])

AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_OUTPUT
