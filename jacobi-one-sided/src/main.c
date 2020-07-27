#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "data.h"

int main(int argc, char** argv)
{
	int rank_world, nprocs_world;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_world);

	instance_t instance;

	if (rank_world == 0)
		read_input(stdin, &instance);

	broadcast_input_data(&instance);


	MPI_Finalize();
	return 0;
}