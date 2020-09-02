# Designed for Intel MPI

MPICC=mpiicc
OPT_FLAGS=-Ofast -fno-alias
STD=-std=c99
OUTNAME=jacobi2d.x

all: release

debug:
	$(MPICC) -DDEBUG -g -O0 -Wall -Iinclude -o $(OUTNAME) src/*.c

release:
	$(MPICC) $(STD) $(OPT_FLAGS) -xHost -Iinclude -o $(OUTNAME) src/*.c

ivybridge:
	$(MPICC) $(STD) $(OPT_FLAGS) -xIVYBRIDGE -Iinclude -o $(OUTNAME) src/*.c

skylake:
	$(MPICC) $(STD) $(OPT_FLAGS) -xSKYLAKE -Iinclude -o $(OUTNAME) src/*.c
