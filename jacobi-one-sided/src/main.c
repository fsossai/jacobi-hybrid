#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "data.h"
#include "kernel.h"
#include "arguments.h"

int main(int argc, char* argv[])
{
	int rank_world, nprocs_world;

	instance_t instance;
	memset(&instance, 0x00, sizeof(instance_t));

	if (parse_command_line_arguments(argc, argv, &instance) == ERROR_PARSING_ARGUMENTS)
		return -1;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_world);

	/*if (rank_world == 0)
		read_input(stdin, &instance);*/
	if (rank_world == 0) // debug
	{
		instance.domain_sizes[0] = 300;
		instance.domain_sizes[1] = 300;
		instance.domain_sizes[2] = 300;
		instance.alpha = 0.8;
		instance.relaxation = 1.0;
		instance.tolerance = 1e-10;
		instance.max_iterations = 100;
	}

	// printing input data
	if (rank_world == 0)
	{
		printf("Domain size\t\t: %ix%ix%i\n",
			instance.domain_sizes[0],
			instance.domain_sizes[1],
			instance.domain_sizes[2]);
		printf("Alpha\t\t\t: %.5lf\n", instance.alpha);
		printf("Relaxation\t\t: %.5lf\n", instance.relaxation);
		printf("Tolerance\t\t: %e\n", instance.tolerance);
		printf("Max iteration\t\t: %i\n", instance.max_iterations);
	}

	// creating shared and head communicators
	MPI_Comm comm_shared, comm_head;
	setup_shared_and_heads(&instance, &comm_shared, &comm_head);
	int nprocs_head;
	if (comm_head != MPI_COMM_NULL)
		MPI_Comm_size(comm_head, &nprocs_head);

	// to keep things simple, min{dims_i} must be >= 'nprocs'
	broadcast_input_data_head(comm_head, &instance);

	// creating a cartesian topology upon 'comm_head'
	MPI_Comm comm_cart = MPI_COMM_NULL;
	int coords[DOMAIN_DIM];
	setup_topology(comm_head, instance.cart_splits, coords, &comm_cart);

	// computing global and local subdomains' offsets and sizes
	compute_subdomains(comm_cart, coords, instance.cart_splits, &instance);
	broadcast_data_shared(comm_shared, &instance);
	compute_local_workload(comm_shared, &instance);
	initialize_problem(comm_cart, &instance);
	allocate_shared_resources(comm_cart, comm_shared, &instance);

	double local_timer = -MPI_Wtime();
	compute_jacobi(comm_cart, comm_shared, &instance);
	local_timer += MPI_Wtime();

	#ifdef DEBUG
	for (int p = 0; p < nprocs_world; p++)
	{
		if (rank_world == p)
		{
			int rank_shared;
			MPI_Comm_rank(comm_shared, &rank_shared);
			printf("w%2i s%2i c(%2i,%2i,%2i) ",
				rank_world, rank_shared,
				coords[0], coords[1], coords[2]);
			printf("sdsizes %5i  %5i %5i, offs (%5i+%4i,%5i+%4i,%5i+%4i) lsdsizes (%4i %4i %4i)\n",
				instance.subdomain_sizes[0],
				instance.subdomain_sizes[1],
				instance.subdomain_sizes[2],
				instance.subdomain_offsets[0],
				instance.local_subdomain_offsets[0],
				instance.subdomain_offsets[1],
				instance.local_subdomain_offsets[1],
				instance.subdomain_offsets[2],
				instance.local_subdomain_offsets[2],
				instance.local_subdomain_sizes[0],
				instance.local_subdomain_sizes[1],
				instance.local_subdomain_sizes[2]);
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	#endif

	if (comm_head != MPI_COMM_NULL)
	{
		double iteration_time_avg;
		MPI_Reduce(&instance.total_computation_time, &iteration_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, comm_head);
		iteration_time_avg /= (double)nprocs_head * (double)instance.performed_iterations;

		if (rank_world == 0)
		{
			printf("Performed iterations\t\t\t: %i\n", instance.performed_iterations);
			printf("Residual\t\t\t\t: %e\n", instance.residual);
			printf("---------------------------------------------------------------\n");
			printf("Using shared memory\t\t\t: %s\n", (instance.use_shared_memory ? "Yes" : "No"));
			if (instance.use_shared_memory)
				printf("Shared subdomain split direction\t: %i\n", instance.local_subdomain_split_direction);
			printf("Cartesian topology arrangement\t\t: %ix%ix%i\n",
				instance.cart_splits[0],
				instance.cart_splits[1],
				instance.cart_splits[2]);
			printf("Total elapsed time\t\t\t: %.3lf s\n", local_timer);
			printf("Average time per iteration\t\t: %.3lf ms\n", iteration_time_avg * 1e3);
			printf("Performance\t\t\t\t: %.3lf MFlops\n",
				(double)instance.performed_iterations *
				(double)instance.domain_sizes[0] *
				(double)instance.domain_sizes[1] *
				(double)instance.domain_sizes[2] /
				local_timer * 1e-6 * (double)KERNEL_FLOPS);
		}
	}

	close_problem(&instance);
	MPI_Finalize();
	return 0;
}