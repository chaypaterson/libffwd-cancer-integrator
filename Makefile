CC = g++
SRC = five-stage-characteristic.cpp
BUILDDIR = bin
BIN = five-step-characteristic
FLAGS = -lm
OPT1 = -O3
OPT2 = $(OPT1) -march=znver3 -flto

fivestage :
	mkdir -p $(BUILDDIR)
	$(CC) $(SRC) $(FLAGS) -o $(BUILDDIR)/$(BIN)
	$(CC) $(SRC) $(FLAGS) $(OPT1) -o $(BUILDDIR)/$(BIN)-opt1
	$(CC) $(SRC) $(FLAGS) $(OPT2) -o $(BUILDDIR)/$(BIN)-opt2
