#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "data.h"
#include "kernel.h"
#include "arguments.h"
#define MASTER 0

void print_input_info(instance_t* instance);
void print_debug_info(instance_t* instance, int* coords, MPI_Comm comm_shared);
void print_configuration(instance_t* instance, MPI_Comm comm_head);
void print_stats(instance_t* instance, MPI_Comm comm_head);

int main(int argc, char* argv[])
{
	int rank_world, nprocs_world;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_world);

	instance_t instance;
	memset(&instance, 0x00, sizeof(instance_t));

	int args_error = 0;
	if (rank_world == MASTER)
		args_error = parse_command_line_arguments(argc, argv, &instance);

	MPI_Bcast(&args_error, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	if (args_error == ERROR_PARSING_ARGUMENTS)
	{
		MPI_Finalize();
		return (rank_world == MASTER) ? ERROR_PARSING_ARGUMENTS : 0;
	}

	if (rank_world == MASTER)
		read_input(instance.input_stream, &instance);

	// printing input data
	print_input_info(&instance);

	// creating shared and head communicators
	MPI_Comm comm_shared, comm_head;
	setup_shared_and_heads(&instance, &comm_shared, &comm_head);

	// to keep things simple, min{dims_i} must be >= 'nprocs'
	broadcast_input_data_head(comm_head, &instance);

	// creating a cartesian topology upon 'comm_head'
	MPI_Comm comm_cart = MPI_COMM_NULL;
	int coords[DOMAIN_DIM];
	setup_topology(comm_head, instance.cart_splits, coords, &comm_cart);

	// computing global and local subdomains' offsets and sizes
	compute_subdomains(comm_cart, coords, instance.cart_splits, &instance);

	// some input data are broadcasted inside every shared memory islands
	broadcast_data_shared(comm_shared, &instance);

	// calculating how to split the work among the processes in the shared region
	compute_local_workload(comm_shared, &instance);

	// matrices are allocated and F is initialized
	initialize_problem(comm_cart, &instance);
	allocate_shared_resources(comm_cart, comm_shared, &instance);
	print_configuration(&instance, comm_head);

	#ifdef DEBUG
	print_debug_info(&instance, coords, comm_shared);
	#endif

	compute_jacobi(comm_cart, comm_shared, &instance);

	print_stats(&instance, comm_head);

	close_problem(&instance);
	MPI_Finalize();
	return 0;
}

void print_input_info(instance_t* instance)
{
	int rank_world;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);

	if (rank_world == MASTER)
	{
		printf("Domain size\t\t: %ix%ix%i\n",
			instance->domain_sizes[0],
			instance->domain_sizes[1],
			instance->domain_sizes[2]);
		printf("Alpha\t\t\t: %.5lf\n", instance->alpha);
		printf("Relaxation\t\t: %.5lf\n", instance->relaxation);
		printf("Tolerance\t\t: %e\n", instance->tolerance);
		printf("Max iteration\t\t: %i\n", instance->max_iterations);
		printf("---------------------------------------------------------------\n");
	}
}

void print_debug_info(instance_t *instance, int *coords, MPI_Comm comm_shared)
{
	int rank_shared, rank_world, nprocs_world;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_world);
	MPI_Comm_rank(comm_shared, &rank_shared);

	for (int p = 0; p < nprocs_world; p++)
	{
		if (rank_world == p)
		{
			printf("w%2i s%2i c(%2i,%2i,%2i) ",
				rank_world, rank_shared,
				coords[0], coords[1], coords[2]);
			printf("sdsizes (%5i,%5i,%5i), offs (%5i+%4i,%5i+%4i,%5i+%4i) lsdsizes (%4i,%4i,%4i)\n",
				instance->subdomain_sizes[0],
				instance->subdomain_sizes[1],
				instance->subdomain_sizes[2],
				instance->subdomain_offsets[0],
				instance->local_subdomain_offsets[0],
				instance->subdomain_offsets[1],
				instance->local_subdomain_offsets[1],
				instance->subdomain_offsets[2],
				instance->local_subdomain_offsets[2],
				instance->local_subdomain_sizes[0],
				instance->local_subdomain_sizes[1],
				instance->local_subdomain_sizes[2]);
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void print_configuration(instance_t* instance, MPI_Comm comm_head)
{
	int rank_world, nprocs_world;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_world);

	if (rank_world == MASTER)
	{
		printf("Total number of processes\t\t: %i\n", nprocs_world);
		printf("Using shared memory\t\t\t: %s\n", (instance->use_shared_memory ? "Yes" : "No"));
		if (instance->use_shared_memory)
			printf("Shared subdomain split direction\t: %i\n", instance->local_subdomain_split_direction);
		printf("Cartesian topology arrangement\t\t: %ix%ix%i\n",
			instance->cart_splits[0],
			instance->cart_splits[1],
			instance->cart_splits[2]);
		printf("Num heads per physical shared region\t: %i\n", instance->heads_per_shared_region);
	}
}

void print_stats(instance_t* instance, MPI_Comm comm_head)
{
	int rank_world, nprocs_head;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_head);

	if (comm_head != MPI_COMM_NULL)
	{
		double iteration_time_avg;
		MPI_Reduce(&instance->total_computation_time, &iteration_time_avg, 1, MPI_DOUBLE, MPI_SUM, 0, comm_head);
		iteration_time_avg /= (double)nprocs_head * (double)instance->performed_iterations;

		if (rank_world == MASTER)
		{
			printf("---------------------------------------------------------------\n");
			printf("Performed iterations\t\t\t: %i\n", instance->performed_iterations);
			printf("Residual\t\t\t\t: %e\n", instance->residual);
			printf("Total elapsed time\t\t\t: %.3lf s\n", instance->total_computation_time);
			printf("Average time per iteration\t\t: %.3lf ms\n", iteration_time_avg * 1e3);
			printf("Performance\t\t\t\t: %.3lf MFlops\n",
				(double)instance->performed_iterations *
				(double)instance->domain_sizes[0] *
				(double)instance->domain_sizes[1] *
				(double)instance->domain_sizes[2] /
				instance->total_computation_time * 1e-6 * (double)KERNEL_FLOPS);
		}
	}
}