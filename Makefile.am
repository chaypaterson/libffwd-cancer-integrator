AUTOMAKE_OPTIONS = subdir-objects foreign
dist_doc_DATA = README.md

STD = --std=c++17
OPT1 = -O3

AM_CPPFLAGS = -I$(src_topdir)../include

## Core libraries:
lib_LTLIBRARIES = libffwd.la
libffwd_la_SOURCES = src/core/fast-forward.cpp
libffwd_la_LDFLAGS = -version-info 0:0:0

## Unit tests:
bin_PROGRAMS = unittest
unittest_SOURCES = src/tests/likelihood-unit-test.cpp
unittest_LDADD = libffwd.la
unittest_CPPFLAGS = $(AM_CPPFLAGS)

## Error estimation:

## Inference harness:

## Quality of life:
clean-local:
	@rm config.status configure config.log
	@rm Makefile
	@rm -r autom4te.cache/
	@rm aclocal.m4
	@rm compile install-sh missing Makefile.in
