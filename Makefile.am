dist_doc_DATA = README.md

STD = --std=c++17
OPT1 = -O3

AM_CPPFLAGS = -I$(top_srcdir)/core/

bin_PROGRAMS = gillespie_unit_test likelihood_unit_test
gillespie_unit_test_SOURCES = src/tests/gillespie-sampler.cpp
gillespie_unit_test_CPPFLAGS = -Isrc/core
gillespie_unit_test_CXXFLAGS = -lgsl -lgslcblas -pthread $STD $OPT1
gillespie_unit_test_LDADD = libgillespie.la

likelihood_unit_test_SOURCES = src/tests/likelihood-unit-test.cpp
likelihood_unit_test_LDADD = libffwd.la

lib_LTLIBRARIES = libffwd.la libgillespie.la
include_HEADERS = src/core/fast-forward.hpp src/core/graph-model-spec.hpp src/core/gillespie-algorithm.hpp

libffwd_la_SOURCES = src/core/fast-forward.cpp src/core/fast-forward.hpp src/core/graph-model-spec.hpp
libffwd_la_CPPFLAGS = -Isrc/core
libffwd_la_CXXFLAGS = $STD $OPT1

libgillespie_la_SOURCES = src/core/gillespie-algorithm.cpp
libgillespie_la_LIBADD = -lgsl -lgslcblas

