AUTOMAKE_OPTIONS = subdir-objects foreign
dist_doc_DATA = README.md

AM_CPPFLAGS = -I$(src_topdir)../include

## Core libraries:
lib_LTLIBRARIES = libffwd.la libgillespie.la
libffwd_la_SOURCES = src/core/fast-forward.cpp
libffwd_la_LDFLAGS = -version-info 0:0:0
libgillespie_la_SOURCES = src/core/gillespie-algorithm.cpp
libgillespie_la_LDFLAGS = -version-info 0:0:0

## Unit tests:
tests:
	test -d tests || mkdir tests
bin_testsdir = $(bindir)/tests
bin_tests_PROGRAMS = likelihoodunittest gillespiesampler \
                        fivestage tslossffwd tslossgill
likelihoodunittest_SOURCES = src/tests/likelihood-unit-test.cpp
likelihoodunittest_LDADD = libffwd.la
likelihoodunittest_CPPFLAGS = $(AM_CPPFLAGS)
gillespiesampler_SOURCES = src/tests/gillespie-sampler.cpp
gillespiesampler_LDADD = libgillespie.la
fivestage_SOURCES = src/tests/five-stage-characteristic.cpp
fivestage_LDADD = libffwd.la
tslossffwd_SOURCES = src/tests/ts-loss-characteristic.cpp
tslossffwd_LDADD = libffwd.la
tslossgill_SOURCES = src/tests/ts-loss-gillespie.cpp
tslossgill_LDADD = libgillespie.la

bin_PROGRAMS = tslosserrs gilllosserrs guesser
## Error estimation:
tslosserrs_SOURCES = src/errors/ts-loss-characteristic-errors.cpp
tslosserrs_LDADD = libffwd.la
gilllosserrs_SOURCES = src/errors/gillespie-errors.cpp
gilllosserrs_LDADD = libgillespie.la
## Inference harness:
guesser_SOURCES = src/inference/likelihood-optimisation.cpp
guesser_LDADD = libffwd.la libgillespie.la
guesser_CPPFLAGS = $(AM_CPPFLAGS)

## Clean the build directory:
clean-local:
	@rm Makefile
	@rm -r autom4te.cache/
	@rm aclocal.m4
	@rm compile install-sh missing Makefile.in
