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
	const int nheads_per_node = 1;
	setup_shared_and_heads(nheads_per_node, &comm_shared, &comm_head);
	int nprocs_head;
	if (comm_head != MPI_COMM_NULL)
		MPI_Comm_size(comm_head, &nprocs_head);

	// to keep things simple, min{dims_i} must be >= 'nprocs'
	broadcast_input_data_head(comm_head, &instance);

	// creating a cartesian topology upon 'comm_head'
	MPI_Comm comm_cart = MPI_COMM_NULL;
	int nsplits_per_dim[DOMAIN_DIM], coords[DOMAIN_DIM];
	setup_topology(comm_head, nsplits_per_dim, coords, &comm_cart);

	// computing global and local subdomains' offsets and sizes
	compute_subdomains(comm_cart, coords, nsplits_per_dim, &instance);
	broadcast_data_shared(comm_shared, &instance);
	compute_local_workload(comm_shared, &instance);
	initialize_problem(comm_cart, &instance);
	allocate_shared_resources(comm_cart, comm_shared, &instance);

	double local_timer = -MPI_Wtime();
	compute_jacobi(comm_cart, comm_shared, &instance);
	local_timer += MPI_Wtime();

	const char show = 1;
	for (int p = 0; show && p < nprocs_world; p++)
	{
		if (rank_world == p)
		{
			int rank_shared = 666; //should never be displayed
			MPI_Comm_rank(comm_shared, &rank_shared);
			printf("w%2i s%2i c(%2i,%2i,%2i) ",
				rank_world, rank_shared,
				coords[0], coords[1], coords[2]);
			printf("sd sizes %5i  %5i %5i, offs (%5i+%4i,%5i+%4i,%5i+%4i) lsdsizes (%4i %4i %4i)\n",
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
			//printf(" alpha: %5.2lf, maxit %i, tol %lf relax %lf\n",
			//	instance.alpha, instance.max_iterations, instance.tolerance, instance.relaxation);
			//if (comm_head != MPI_COMM_NULL)
			//	print_subdomain(instance.U, &instance, "%8.3lf ");
			//print_F(&instance, "%8.3lf ");
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (comm_head != MPI_COMM_NULL)
	{
		double iteration_time_avg;
		MPI_Reduce(&instance.total_computation_time, &iteration_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, comm_head);
		iteration_time_avg /= (double)nprocs_head * (double)instance.performed_iterations;

		if (rank_world == 0)
		{
			int shared_mem = 1;
			#ifdef NO_SHARED_MEMORY
			shared_mem = 0;
			#endif 

			printf("Performed iterations\t\t\t: %i\n", instance.performed_iterations);
			printf("Residual\t\t\t\t: %e\n", instance.residual);
			printf("---------------------------------------------------------------\n");
			printf("Using shared memory\t\t\t: %s\n", (shared_mem ? "Yes" : "No"));
			if (shared_mem)
				printf("Shared subdomain split direction\t: %i\n", instance.local_subdomain_split_direction);
			printf("Cartesian topology arrangement\t\t: %ix%ix%i\n",
				nsplits_per_dim[0],
				nsplits_per_dim[1],
				nsplits_per_dim[2]);
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