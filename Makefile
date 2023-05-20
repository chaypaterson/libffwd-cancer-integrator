CC = g++
CHARSRC = five-stage-characteristic.cpp
BUILDDIR = bin
LIBDIR = libs
BIN = five-step-characteristic
OPT1 = -O3
STD = --std=c++11
MACLIBS = -I/opt/homebrew/include/ -L/opt/homebrew/lib/
#OPTRYZEN = $(OPT1) -march=znver3 -flto
LIBCHAR = $(LIBDIR)/libflying.so
LIBGILL = $(LIBDIR)/libgillespie.so
FLAGS = $(STD) $(MACLIBS)
GILLFLAGS = $(MACLIBS) -lgsl -lgslcblas -pthread

BIN2HIT = two-hit-characteristic

#TODO better names for build dir and src dirs. move core of source to src/core/.

guesser : libflying.so libgillespie.so builddir
	$(CC) $(LIBCHAR) $(LIBGILL) parameter_inference/likelihood-optimisation.cpp $(STD) $(GILLFLAGS) $(OPT1) -o $(BUILDDIR)/guesser

sampler : libgillespie.so builddir
	$(CC) $(LIBGILL) parameter_inference/gillespie-sampler.cpp $(STD) $(GILLFLAGS) $(OPT1) -o $(BUILDDIR)/gillespie_sampler

fivestage : $(CHARSRC) libflying.so builddir
	$(CC) $(LIBCHAR) src/tests/$(CHARSRC) $(FLAGS) -o $(BUILDDIR)/$(BIN)
	$(CC) $(LIBCHAR) src/tests/$(CHARSRC) $(FLAGS) $(OPT1) -o $(BUILDDIR)/$(BIN)-opt1

twostage : $(CHARSRC) libflying.so libgillespie.so builddir
	$(CC) $(LIBCHAR) $(BIN2HIT).cpp $(FLAGS) $(OPT2) -o $(BUILDDIR)/$(BIN2HIT)
	$(CC) $(LIBGILL) two-hit-gillespie.cpp $(GILLFLAGS) $(OPT2) -o $(BUILDDIR)/two-hit-gillespie-2

tsloss : libflying.so libgillespie.so builddir
	$(CC) $(LIBCHAR) ts-loss-characteristic.cpp $(FLAGS) $(OPT1) -o $(BUILDDIR)/tsconj
	$(CC) $(LIBGILL) ts-loss-gillespie.cpp $(STD) $(GILLFLAGS) $(OPT1) -o $(BUILDDIR)/tsgillespie

unittests : libflying.so libgillespie.so builddir
	$(CC) $(LIBGILL) src/tests/gillespie-sampler.cpp $(GILLFLAGS) -std=c++11 $(OPT1) -o bin/gillespie_sampler
	$(CC) $(LIBCHAR) src/tests/likelihood-unit-test.cpp $(FLAGS) $(OPT1) -o bin/unittest

libflying.so: libs
	$(CC) flying-conjugates.cpp $(FLAGS) $(OPT1) -c -o $(LIBCHAR)

libgillespie.so: libs
	$(CC) gillespie-algorithm.cpp $(FLAGS) $(OPT1) -c -o $(LIBGILL)

libs :
	mkdir -p $(LIBDIR)

builddir :
	mkdir -p $(BUILDDIR)

