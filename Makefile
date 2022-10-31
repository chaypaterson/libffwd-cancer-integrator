CC = g++
CHARSRC = five-stage-characteristic.cpp
BUILDDIR = bin
LIBDIR = libs
BIN = five-step-characteristic
FLAGS = -lm
OPT1 = -O3
#OPT2 = $(OPT1) -march=znver3 -flto
LIBCHAR = $(LIBDIR)/libflying.so
LIBGILL = $(LIBDIR)/libgillespie.so
GILLFLAGS = -lgsl -pthread

BIN2HIT = two-hit-characteristic

libs :
	mkdir -p $(LIBDIR)

libflying.so: libs
	$(CC) flying-conjugates.cpp $(FLAGS) $(OPT1) -c -o $(LIBCHAR)

libgillespie.so: libs
	$(CC) gillespie-algorithm.cpp $(FLAGS) $(OPT1) -c -o $(LIBGILL)

fivestage : $(CHARSRC) libflying.so
	mkdir -p $(BUILDDIR)
	$(CC) $(LIBCHAR) $(CHARSRC) $(FLAGS) -o $(BUILDDIR)/$(BIN)
	$(CC) $(LIBCHAR) $(CHARSRC) $(FLAGS) $(OPT1) -o $(BUILDDIR)/$(BIN)-opt1
	$(CC) $(LIBCHAR) $(CHARSRC) $(FLAGS) $(OPT2) -o $(BUILDDIR)/$(BIN)-opt2

twostage : $(CHARSRC) libflying.so libgillespie.so
	mkdir -p $(BUILDDIR)
	$(CC) $(LIBCHAR) $(BIN2HIT).cpp $(FLAGS) $(OPT2) -o $(BUILDDIR)/$(BIN2HIT)
	$(CC) $(LIBGILL) two-hit-gillespie.cpp $(GILLFLAGS) $(OPT2) -o $(BUILDDIR)/two-hit-gillespie-2

tsloss :
	$(CC) $(LIBCHAR) ts-loss-characteristic.cpp $(FLAGS) $(OPT1) -o $(BUILDDIR)/tsconj
	$(CC) $(LIBGILL) ts-loss-gillespie.cpp $(GILLFLAGS) $(OPT1) -o $(BUILDDIR)/tsgillespie
