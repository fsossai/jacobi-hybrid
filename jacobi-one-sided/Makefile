# Designed for Intel MPI

MPICC = mpiicc

release:
	$(MPICC) -Ofast -Iinclude -ojacobi src/*.c

debug:
	$(MPICC) -g -O0 -Iinclude -ojacobi src/*.c
