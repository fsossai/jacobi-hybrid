#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "data.h"
#include "kernel.h"

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
		instance.domain_sizes[0] = 10;
		instance.domain_sizes[1] = 15;
		instance.alpha = 0.8;
		instance.relaxation = 1.0;
		instance.tolerance = 1e-16;
		instance.max_iterations = 10;
	}

	// creating shared and head communicators
	MPI_Comm comm_shared, comm_head;
	const int nheads_per_node = 1;
	setup_shared_and_heads(nheads_per_node, &comm_shared, &comm_head);

	// to keep things simple, min{dims_i} must be >= 'nprocs'
	broadcast_input_data_head(comm_head, &instance);	

	// creating a cartesian topology upon 'comm_head'
	MPI_Comm comm_cart = MPI_COMM_NULL;
	int nprocs_per_dim[DOMAIN_DIM], coords[DOMAIN_DIM];
	setup_topology(comm_head, nprocs_per_dim, coords, &comm_cart);
	
	// computing subdomain offsets and sizes
	compute_limits(comm_cart, coords, nprocs_per_dim, &instance);
	broadcast_data_shared(comm_shared, &instance);

	initialize_problem(comm_cart, &instance);
	
	compute_jacobi(comm_cart, &instance);

	const char show = 1;
	for (int p = 0; show && p < nprocs_world; p++)
	{
		if (rank_world == p)
		{
			int rank_shared = 666; //should never be displayed
			MPI_Comm_rank(comm_shared, &rank_shared);
			printf("w%2i s%2i c(%2i,%2i) ",
				rank_world, rank_shared,
				coords[0], coords[1]);
			printf("sd sizes %2i %2i, offs %3i %3i\n",
				instance.subdomain_sizes[0],
				instance.subdomain_sizes[1],
				instance.subdomain_offsets[0],
				instance.subdomain_offsets[1]);
			printf(" alpha: %5.2lf, maxit %i, tol %lf relax %lf\n",
				instance.alpha, instance.max_iterations, instance.tolerance, instance.relaxation);
			if (rank_shared == 0)
				print_subdomain(instance.U, &instance, "%8.3lf ");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	close_problem(&instance);
	MPI_Finalize();
	return 0;
}