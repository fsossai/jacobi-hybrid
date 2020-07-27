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
		instance.domain_sizes[0] = 100;
		instance.domain_sizes[1] = 100;
		instance.domain_sizes[2] = 100;
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
	int rank_shared, nprocs_shared;
	int rank_head, nprocs_head, color;
	rank_head = MPI_PROC_NULL; //D
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shared);
	MPI_Comm_rank(comm_shared, &rank_shared);
	MPI_Comm_size(comm_shared, &nprocs_shared);
	color = (rank_shared == 0) ? 1 : MPI_UNDEFINED;
	MPI_Comm_split(comm_shared, color, 0, &comm_head);
	if (comm_head != MPI_COMM_NULL)
	{
		MPI_Comm_rank(comm_head, &rank_head);
		MPI_Comm_size(comm_head, &nprocs_head);
	}

	// creating a cartesian topology upon 'comm_head'
	MPI_Comm comm_cart = MPI_COMM_NULL;
	int rank_cart = MPI_PROC_NULL;
	int nprocs_per_dim[DOMAIN_DIM], periods[DOMAIN_DIM];
	int coords[DOMAIN_DIM];
	coords[0] = coords[1] = coords[2] = -1; //D
	memset(periods, 0x00, DOMAIN_DIM * sizeof(int));
	MPI_Dims_create(nprocs_head, DOMAIN_DIM, nprocs_per_dim);
	if (comm_head != MPI_COMM_NULL)
		MPI_Cart_create(comm_head, DOMAIN_DIM, nprocs_per_dim, periods, 0, &comm_cart);
	if (comm_cart != MPI_COMM_NULL)
	{
		MPI_Comm_rank(comm_cart, &rank_cart);
		MPI_Cart_coords(comm_cart, rank_cart, DOMAIN_DIM, coords);
	}

	// computing offsets, and subdomains sizes
	if (comm_cart != MPI_COMM_NULL)
	{
		int divisor, reminder;
		for (int i = 0, index; i < DOMAIN_DIM; i++)
		{
			divisor = instance.domain_sizes[i] / nprocs_per_dim[i];
			reminder = instance.domain_sizes[i] % nprocs_per_dim[i];
			index = coords[i];
			instance.subdomain_offsets[i] = index * divisor + ((index <= reminder) ? index : reminder);
			index = coords[i] + 1;
			instance.subdomain_sizes[i] =
				index * divisor + ((index <= reminder) ? index : reminder) -
				instance.subdomain_offsets[i];
		}
		initialize_problem(&instance);
	}

	for (int p = 0; p < nprocs_world; p++)
	{
		if (rank_world == p)
		{
			printf("w%2i s%2i h%2i c(%2i,%2i,%2i) ",
				rank_world, rank_shared, rank_head,
				coords[0], coords[1], coords[2]);
			printf("sd sizes %2i %2i %2i, offs %2i %2i %2i\n",
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