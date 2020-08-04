# Designed for Intel MPI

MPICC = mpiicc

all: clean release

debug:
	$(MPICC) -g -O0 -Wall -Iinclude -ojacobi.x src/*.c

release:
	$(MPICC) -Ofast -xAVX -fno-alias -Iinclude -ojacobi.x src/*.c

clean:
	rm -f jacobi.x