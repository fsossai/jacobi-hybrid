# Designed for Intel MPI

MPICC = mpiicc

debug:
	$(MPICC) -g -O0 -Wall -Iinclude -ojacobi.x src/*.c

all: clean release

release:
	$(MPICC) -Ofast -xAVX -fno-alias -Iinclude -ojacobi.x src/*.c

clean:
	rm -f *.x