AUTOMAKE_OPTIONS = subdir-objects foreign
dist_doc_DATA = README.md

AM_CPPFLAGS = -I$(src_topdir)../include

## Core libraries:
lib_LTLIBRARIES = libffwd.la libgillespie.la
libffwd_la_SOURCES = src/core/fast-forward.cpp \
                     include/graph-model-spec.hpp \
                     include/fast-forward.hpp
libffwd_la_LDFLAGS = -version-info 0:0:0
libgillespie_la_SOURCES = src/core/gillespie-algorithm.cpp \
                          include/graph-model-spec.hpp \
                          include/gillespie-algorithm.hpp
libgillespie_la_LDFLAGS = $(GSL_LIBS) -version-info 0:0:0

## Unit tests:
tests:
	test -d tests || mkdir tests
bin_testsdir = $(bindir)/tests
bin_tests_PROGRAMS = likelihoodunittest gillespiesampler \
                        fivestage tslossffwd tslossgill \
                        tslossgillpoly
likelihoodunittest_SOURCES = src/tests/likelihood-unit-test.cpp \
                             src/tests/max-likelihood.hpp
likelihoodunittest_LDADD = libffwd.la
likelihoodunittest_CPPFLAGS = $(AM_CPPFLAGS)
gillespiesampler_SOURCES = src/tests/gillespie-sampler.cpp
gillespiesampler_LDADD = libgillespie.la
gillespiesampler_LDFLAGS = $(GSL_LIBS)
fivestage_SOURCES = src/tests/five-stage-characteristic.cpp
fivestage_LDADD = libffwd.la
tslossffwd_SOURCES = src/tests/ts-loss-characteristic.cpp
tslossffwd_LDADD = libffwd.la
tslossgill_SOURCES = src/tests/ts-loss-gillespie.cpp
tslossgill_LDADD = libgillespie.la
tslossgill_LDFLAGS = $(GSL_LIBS)
tslossgillpoly_SOURCES = src/tests/ts-loss-gillespie-poly.cpp
tslossgillpoly_LDADD = libgillespie.la
tslossgillpoly_LDFLAGS = $(GSL_LIBS)

bin_PROGRAMS = tslosserrs gilllosserrs guesser
## Error estimation:
tslosserrs_SOURCES = src/errors/ts-loss-characteristic-errors.cpp
tslosserrs_LDADD = libffwd.la
gilllosserrs_SOURCES = src/errors/gillespie-errors.cpp
gilllosserrs_LDADD = libgillespie.la
gilllosserrs_LDFLAGS = $(GSL_LIBS)
## Inference harness:
guesser_SOURCES = src/inference/likelihood-optimisation.cpp \
                  src/inference/likelihood-optimisation.hpp
guesser_LDADD = libffwd.la libgillespie.la
guesser_CPPFLAGS = $(AM_CPPFLAGS)
