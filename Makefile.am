AUTOMAKE_OPTIONS = subdir-objects
dist_doc_DATA = README.md
lib_LIBRARIES = libffwd.a libgillespie.a

libffwd_a_SOURCES = src/core/fast-forward.cpp
libffwd_a_CPPFLAGS = -Isrc/core
libffwd_a_CXXFLAGS = --std=c++11 -O3

libgillespie_a_SOURCES = src/core/gillespie-algorithm.cpp
libgillespie_a_CPPFLAGS = -Isrc/core
libgillespie_a_CXXFLAGS = --std=c++11 -O3
libgillespie_a_LIBADD = -lgsl -lgslcblas

