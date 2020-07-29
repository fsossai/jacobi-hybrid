# Designed for Intel MPI

MPICC = mpiicc

release:
	$(MPICC) -Ofast -Iinclude -ojacobi.x src/*.c

debug:
	$(MPICC) -g -O0 -Wall -Iinclude -ojacobi.x src/*.c
