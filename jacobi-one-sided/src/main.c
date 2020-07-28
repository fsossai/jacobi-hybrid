#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "data.h"

int main(int argc, char* argv[])
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
		instance.domain_sizes[0] = 13;
		instance.domain_sizes[1] = 14;
		instance.domain_sizes[2] = 15;
		instance.alpha = 0.8;
		instance.relaxation = 1.0;
		instance.tolerance = 1e-16;
	}

	// to keep things simple, min{dims_i} must be >= 'nprocs'
	broadcast_input_data(&instance);

	// debug
	/*instance.domain_sizes[0] = 3;
	instance.domain_sizes[1] = 4;
	instance.domain_sizes[2] = 4;
	instance.subdomain_offsets[0] = 0;
	instance.subdomain_offsets[1] = 0;
	instance.subdomain_offsets[2] = rank_world;
	instance.subdomain_sizes[0] = 3;
	instance.subdomain_sizes[1] = 4;
	instance.subdomain_sizes[2] = 1;*/

	// creating shared and head communicators
	MPI_Comm comm_shared, comm_head;
	int nheads_per_node = 2;
	setup_shared_and_heads(nheads_per_node, &comm_shared, &comm_head);

	// creating a cartesian topology upon 'comm_head'
	MPI_Comm comm_cart = MPI_COMM_NULL;
	int nprocs_per_dim[DOMAIN_DIM], coords[DOMAIN_DIM];
	setup_topology(comm_head, nprocs_per_dim, coords, &comm_cart);
	
	// computing subdomain offsets and sizes
	compute_limits(comm_cart, coords, nprocs_per_dim, &instance);

	initialize_problem(comm_cart, &instance);

	for (int p = 0; p < nprocs_world; p++)
	{
		if (rank_world == p)
		{
			int rank_shared = 666; //should never be displayed
			MPI_Comm_rank(comm_shared, &rank_shared);
			printf("w%2i s%2i c(%2i,%2i,%2i) ",
				rank_world, rank_shared,
				coords[0], coords[1], coords[2]);
			printf("sd sizes %2i %2i %2i, offs %3i %3i %3i\n",
				instance.subdomain_sizes[0],
				instance.subdomain_sizes[1],
				instance.subdomain_sizes[2],
				instance.subdomain_offsets[0],
				instance.subdomain_offsets[1],
				instance.subdomain_offsets[2]);
			//print_subdomain(instance.F, &instance, "%8.3lf ");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	close_problem(&instance);
	MPI_Finalize();
	return 0;
}