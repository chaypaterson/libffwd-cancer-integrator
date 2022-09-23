CC = g++
SRC = five-stage-characteristic.cpp
BUILDDIR = bin
BIN = five-step-characteristic
FLAGS = -lm
OPT1 = -O3
OPT2 = $(OPT1) -march=znver3 -flto
LIB = libflying.so

libflying.so:
	$(CC) flying-conjugates.cpp $(FLAGS) $(OPT1) -c -o $(LIB)

fivestage : $(SRC) libflying.so
	mkdir -p $(BUILDDIR)
	$(CC) $(LIB) $(SRC) $(FLAGS) -o $(BUILDDIR)/$(BIN)
	$(CC) $(LIB) $(SRC) $(FLAGS) $(OPT1) -o $(BUILDDIR)/$(BIN)-opt1
	$(CC) $(LIB) $(SRC) $(FLAGS) $(OPT2) -o $(BUILDDIR)/$(BIN)-opt2
