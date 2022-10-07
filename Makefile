CC = g++
CHARSRC = five-stage-characteristic.cpp
BUILDDIR = bin
LIBDIR = libs
BIN = five-step-characteristic
FLAGS = -lm
OPT1 = -O3
OPT2 = $(OPT1) -march=znver3 -flto
LIBCHAR = $(LIBDIR)/libflying.so
LIBGILL = $(LIBDIR)/libgillespie.so

BIN2HIT = two-hit-characteristic

libflying.so:
	mkdir -p $(LIBDIR)
	$(CC) flying-conjugates.cpp $(FLAGS) $(OPT1) -c -o $(LIBCHAR)

libgillespie.so:
	mkdir -p $(LIBDIR)
	$(CC) gillespie-algorithm.cpp $(FLAGS) $(OPT1) -c -o $(LIBGILL)

fivestage : $(CHARSRC) libflying.so
	mkdir -p $(BUILDDIR)
	$(CC) $(LIBCHAR) $(CHARSRC) $(FLAGS) -o $(BUILDDIR)/$(BIN)
	$(CC) $(LIBCHAR) $(CHARSRC) $(FLAGS) $(OPT1) -o $(BUILDDIR)/$(BIN)-opt1
	$(CC) $(LIBCHAR) $(CHARSRC) $(FLAGS) $(OPT2) -o $(BUILDDIR)/$(BIN)-opt2

twostage : $(CHARSRC) libflying.so
	mkdir -p $(BUILDDIR)
	$(CC) $(LIBCHAR) $(BIN2HIT).cpp $(FLAGS) $(OPT2) -o $(BUILDDIR)/$(BIN2HIT)
