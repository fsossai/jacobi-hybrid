# Designed for Intel MPI

MPICC=mpiicc
OPT_FLAGS=-Ofast -fno-alias

all: release

debug:
	$(MPICC) -DDEBUG -g -O0 -Wall -Iinclude -ojacobi.x src/*.c

release:
	$(MPICC) $(OPT_FLAGS) -xHost -Iinclude -ojacobi.x src/*.c

ivybridge:
	$(MPICC) $(OPT_FLAGS) -xIVYBRIDGE -Iinclude -ojacobi.x src/*.c

skylake:
	$(MPICC) $(OPT_FLAGS) -xSKYLAKE -Iinclude -ojacobi.x src/*.c

