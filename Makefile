CC = g++
BUILDDIR = bin
LIBDIR = libs
CORE = src/core
TESTS = src/tests
STD = --std=c++11
OPT1 = -O3
INCLUDE = -I$(CORE)
FLAGS = $(STD) $(OPT1) $(INCLUDE)
LIBCHAR = $(LIBDIR)/libffwd.so
LIBGILL = $(LIBDIR)/libgillespie.so
GILLFLAGS = $(FLAGS) -lgsl -lgslcblas -pthread

OS = $(shell uname)
ifeq ($(OS),Darwin)
	# Mac only:
	MACLIBS = -I/opt/homebrew/include/ -L/opt/homebrew/lib/
	FLAGS += $(MACLIBS)
	GILLFLAGS += $(MACLIBS)
endif
ifeq ($(OS),Linux)
	# Linux+Ryzen only:
	#OPTRYZEN = $(OPT1) -march=znver3 -flto
endif

all : unittests tsloss fivestage sampler guesser numericalerrors

guesser : libffwd.so libgillespie.so builddir
	$(CC) $(LIBCHAR) $(LIBGILL) src/inference/likelihood-optimisation.cpp $(GILLFLAGS) -o $(BUILDDIR)/guesser

numericalerrors : libffwd.so libgillespie.so builddir
	$(CC) $(LIBCHAR) src/errors/ts-loss-characteristic-errors.cpp $(FLAGS) -o bin/tslosserrs

sampler : libgillespie.so builddir
	$(CC) $(LIBGILL) src/inference/gillespie-sampler.cpp $(GILLFLAGS) -o $(BUILDDIR)/gillespie_sampler

fivestage : libffwd.so builddir
	$(CC) $(LIBCHAR) $(TESTS)/five-stage-characteristic.cpp $(FLAGS) -o $(BUILDDIR)/five-step-characteristic

tsloss : libffwd.so libgillespie.so builddir
	$(CC) $(LIBCHAR) $(TESTS)/ts-loss-characteristic.cpp $(FLAGS) -o $(BUILDDIR)/tsconj
	$(CC) $(LIBGILL) $(TESTS)/ts-loss-gillespie.cpp $(GILLFLAGS) -o $(BUILDDIR)/tsgillespie

unittests : libffwd.so libgillespie.so builddir
	$(CC) $(LIBGILL) $(TESTS)/gillespie-sampler.cpp $(GILLFLAGS) -o bin/gillespie_sampler_test
	$(CC) $(LIBCHAR) $(TESTS)/likelihood-unit-test.cpp $(FLAGS) -o bin/unittest

libffwd.so: libs
	$(CC) $(CORE)/fast-forward.cpp $(FLAGS) -c -o $(LIBCHAR)

libgillespie.so: libs
	$(CC) $(CORE)/gillespie-algorithm.cpp $(FLAGS) -c -o $(LIBGILL)

libs :
	mkdir -p $(LIBDIR)

builddir :
	mkdir -p $(BUILDDIR)
