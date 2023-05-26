CC = g++
CHARSRC = five-stage-characteristic.cpp
BUILDDIR = bin
LIBDIR = libs
BIN = five-step-characteristic
OPT1 = -O3
STD = --std=c++11
FLAGS = $(STD) $(OPT1)
LIBCHAR = $(LIBDIR)/libflying.so
LIBGILL = $(LIBDIR)/libgillespie.so
GILLFLAGS = $(FLAGS) -lgsl -lgslcblas -pthread
ifeq ($(UNAME_S),Darwin)
	# Mac only:
	MACLIBS = -I/opt/homebrew/include/ -L/opt/homebrew/lib/
	FLAGS = $(FLAGS) $(MACLIBS)
	GILLFLAGS = $(MACLIBS) $(GILLFLAGS)
endif
ifeq ($(UNAME_S),Linux)
	# Linux+Ryzen only:
	#OPTRYZEN = $(OPT1) -march=znver3 -flto
endif

#TODO better names for build dir and src dirs. move core of source to src/core/.

guesser : libflying.so libgillespie.so builddir
	$(CC) $(LIBCHAR) $(LIBGILL) parameter_inference/likelihood-optimisation.cpp $(GILLFLAGS) -o $(BUILDDIR)/guesser

sampler : libgillespie.so builddir
	$(CC) $(LIBGILL) parameter_inference/gillespie-sampler.cpp $(GILLFLAGS) -o $(BUILDDIR)/gillespie_sampler

fivestage : $(CHARSRC) libflying.so builddir
	$(CC) $(LIBCHAR) src/tests/$(CHARSRC) $(FLAGS) -o $(BUILDDIR)/$(BIN)

tsloss : libflying.so libgillespie.so builddir
	$(CC) $(LIBCHAR) ts-loss-characteristic.cpp $(FLAGS) -o $(BUILDDIR)/tsconj
	$(CC) $(LIBGILL) ts-loss-gillespie.cpp $(GILLFLAGS) -o $(BUILDDIR)/tsgillespie

unittests : libflying.so libgillespie.so builddir
	$(CC) $(LIBGILL) src/tests/gillespie-sampler.cpp $(GILLFLAGS) -o bin/gillespie_sampler
	$(CC) $(LIBCHAR) src/tests/likelihood-unit-test.cpp $(FLAGS) -o bin/unittest

libflying.so: libs
	$(CC) flying-conjugates.cpp $(FLAGS) -c -o $(LIBCHAR)

libgillespie.so: libs
	$(CC) gillespie-algorithm.cpp $(FLAGS) -c -o $(LIBGILL)

libs :
	mkdir -p $(LIBDIR)

builddir :
	mkdir -p $(BUILDDIR)

