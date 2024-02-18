dist_doc_DATA = README.md

lib_LTLIBRARIES = libffwd.la libgillespie.la
include_HEADERS = src/core/fast-forward.hpp src/core/graph-model-spec.hpp src/core/gillespie-algorithm.hpp

libffwd_la_SOURCES = src/core/fast-forward.cpp src/core/fast-forward.hpp src/core/graph-model-spec.hpp
libffwd_la_CPPFLAGS = -Isrc/core
libffwd_la_CXXFLAGS = --std=c++11 -O3

libgillespie_la_SOURCES = src/core/gillespie-algorithm.cpp
libgillespie_la_LIBADD = -lgsl -lgslcblas

bin_PROGRAMS = gillespie_unit_test likelihood_unit_test
gillespie_unit_test_SOURCES = src/tests/gillespie-sampler.cpp
gillespie_unit_test_CPPFLAGS = -Isrc/core
gillespie_unit_test_CXXFLAGS = -lgsl -lgslcblas -pthread
