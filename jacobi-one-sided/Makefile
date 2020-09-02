# Designed for Intel MPI

MPICC=mpiicc
OPT_FLAGS=-Ofast -fno-alias
OUTNAME=jacobi3d.x

all: release

debug:
	$(MPICC) -DDEBUG -g -O0 -Wall -Iinclude -o $(OUTNAME) src/*.c

release:
	$(MPICC) $(OPT_FLAGS) -xHost -Iinclude -o $(OUTNAME) src/*.c

ivybridge:
	$(MPICC) $(OPT_FLAGS) -xIVYBRIDGE -Iinclude -o $(OUTNAME) src/*.c

skylake:
	$(MPICC) $(OPT_FLAGS) -xSKYLAKE -Iinclude -o $(OUTNAME) src/*.c
