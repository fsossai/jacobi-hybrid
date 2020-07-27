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
	memset(&instance, 0x00, sizeof(instance_t));

	/*if (rank_world == 0)
		read_input(stdin, &instance);*/
	if (rank_world == 0) // debug
	{
		instance.domain_x_size = 100;
		instance.domain_y_size = 100;
		instance.domain_z_size = 100;
		instance.alpha = 0.8;
		instance.relaxation = 1.0;
		instance.tolerance = 1e-16;
	}


	broadcast_input_data(&instance);

	// just for debug purposes
	instance.domain_x_size *= instance.subdomain_x_size = 2;
	instance.domain_y_size *= instance.subdomain_y_size = 3;
	instance.domain_z_size *= instance.subdomain_z_size = 4;

	if (rank_world == 0)
	{
		instance.offset_x = 0;
		instance.offset_y = 0;
		instance.offset_z = 0;
	}
	if (rank_world == 1)
	{
		instance.offset_x = 2;
		instance.offset_y = 3;
		instance.offset_z = 4;
	}

	initialize_problem(&instance);

	/*
	for (int p = 0; p < nprocs_world; p++)
	{
		if (rank_world == p)
		{
			printf("[%i]\n", rank_world);
			print_subdomain(instance.F, &instance, "%8.3lf ");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}*/

	close_problem(&instance);
	MPI_Finalize();
	return 0;
}