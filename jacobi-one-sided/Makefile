# Designed for Intel MPI

MPICC = mpiicc

all: clean release

release:
	$(MPICC) -Ofast -xAVX -fno-alias -Iinclude -ojacobi.x src/*.c

debug:
	$(MPICC) -g -O0 -Wall -Iinclude -ojacobi.x src/*.c

clean:
	rm -f *.x