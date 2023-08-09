CC = g++
BUILDDIR = bin
LIBDIR = libs
CORE = src/core
TESTS = src/tests
NUMER = src/errors
LEARN = src/inference
STD = --std=c++11
OPT1 = -O3
INCLUDE = -I$(CORE)
FLAGS = $(STD) $(OPT1) $(INCLUDE)
LIBFFWD = $(LIBDIR)/libffwd.so
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

guesser : $(LIBFFWD) libgillespie.so builddir
	$(CC) $(LIBFFWD) $(LIBGILL) -I /usr/include/eigen3 $(LEARN)/likelihood-optimisation.cpp $(GILLFLAGS) -o $(BUILDDIR)/guesser

sampler : libgillespie.so builddir
	$(CC) $(LIBGILL) $(LEARN)/gillespie-sampler.cpp $(GILLFLAGS) -o $(BUILDDIR)/gillespie_sampler

numericalerrors : $(LIBFFWD) libgillespie.so builddir
	$(CC) $(LIBFFWD) $(NUMER)/ts-loss-characteristic-errors.cpp $(FLAGS) -o bin/tslosserrs
	$(CC) $(LIBGILL) $(NUMER)/gillespie-errors.cpp $(GILLFLAGS) -o bin/gilllosserrs

tests: fivestage tsloss unittests

fivestage : $(LIBFFWD) builddir
	$(CC) $(LIBFFWD) $(TESTS)/five-stage-characteristic.cpp $(FLAGS) -o $(BUILDDIR)/five-step-characteristic

tsloss : $(LIBFFWD) libgillespie.so builddir
	$(CC) $(LIBFFWD) $(TESTS)/ts-loss-characteristic.cpp $(FLAGS) -o $(BUILDDIR)/tsconj
	$(CC) $(LIBGILL) $(TESTS)/ts-loss-gillespie.cpp $(GILLFLAGS) -o $(BUILDDIR)/tsgillespie

unittests : $(LIBFFWD) libgillespie.so builddir
	$(CC) $(LIBGILL) $(TESTS)/gillespie-sampler.cpp $(GILLFLAGS) -o bin/gillespie_sampler_test
	$(CC) $(LIBFFWD) $(TESTS)/likelihood-unit-test.cpp $(FLAGS) -o bin/unittest

$(LIBFFWD): libs
	$(CC) $(CORE)/fast-forward.cpp $(FLAGS) -c -o $(LIBFFWD)

libgillespie.so: libs
	$(CC) $(CORE)/gillespie-algorithm.cpp $(FLAGS) -c -o $(LIBGILL)

core: $(LIBFFWD) libgillespie.so

libs :
	mkdir -p $(LIBDIR)

builddir :
	mkdir -p $(BUILDDIR)
