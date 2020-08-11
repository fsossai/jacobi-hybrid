#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"

void read_input(instance_t* instance)
{
	FILE* input = instance->input_stream;
	if (input == NULL)
	{
		input = stdin;
		printf("Reading instance data from standard input.\n");
	}

	for (int i = 0; i < DOMAIN_DIM; i++)
		fscanf(input, "%i", &instance->domain_sizes[i]);
	fscanf(input, "%lf", &instance->alpha);
	fscanf(input, "%lf", &instance->relaxation);
	fscanf(input, "%lf", &instance->tolerance);
	fscanf(input, "%i", &instance->max_iterations);
	fclose(input);
}

void broadcast_input_data_head(MPI_Comm comm_head, instance_t * instance)
{
	MPI_Aint displacements[7];
	displacements[0] = (MPI_Aint)offsetof(instance_t, domain_sizes);
	displacements[1] = (MPI_Aint)offsetof(instance_t, cart_splits);
	displacements[2] = (MPI_Aint)offsetof(instance_t, local_subdomain_split_direction);
	displacements[3] = (MPI_Aint)offsetof(instance_t, alpha);
	displacements[4] = (MPI_Aint)offsetof(instance_t, relaxation);
	displacements[5] = (MPI_Aint)offsetof(instance_t, tolerance);
	displacements[6] = (MPI_Aint)offsetof(instance_t, max_iterations);

	int block_lengths[7];
	block_lengths[0] = DOMAIN_DIM;
	block_lengths[1] = DOMAIN_DIM;
	for (int i = 2; i < 7; i++)
		block_lengths[i] = 1;

	MPI_Datatype types[7] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };

	MPI_Datatype parameters_struct;
	MPI_Type_create_struct(7, block_lengths, displacements, types, &parameters_struct);
	MPI_Type_commit(&parameters_struct);

	int root;
	if (comm_head != MPI_COMM_NULL)
		MPI_Bcast(instance, 1, parameters_struct, root = 0, comm_head);

	MPI_Type_free(&parameters_struct);
}

void broadcast_data_shared(MPI_Comm comm_shared, instance_t * instance)
{
	MPI_Aint displacements[8];
	displacements[0] = (MPI_Aint)offsetof(instance_t, domain_sizes);
	displacements[1] = (MPI_Aint)offsetof(instance_t, subdomain_sizes);
	displacements[2] = (MPI_Aint)offsetof(instance_t, subdomain_offsets);
	displacements[3] = (MPI_Aint)offsetof(instance_t, local_subdomain_split_direction);
	displacements[4] = (MPI_Aint)offsetof(instance_t, alpha);
	displacements[5] = (MPI_Aint)offsetof(instance_t, relaxation);
	displacements[6] = (MPI_Aint)offsetof(instance_t, tolerance);
	displacements[7] = (MPI_Aint)offsetof(instance_t, max_iterations);

	int block_lengths[8];
	block_lengths[0] = DOMAIN_DIM;
	block_lengths[1] = DOMAIN_DIM;
	block_lengths[2] = DOMAIN_DIM;
	for (int i = 3; i < 8; i++)
		block_lengths[i] = 1;

	MPI_Datatype types[8] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };

	MPI_Datatype parameters_struct;
	MPI_Type_create_struct(8, block_lengths, displacements, types, &parameters_struct);
	MPI_Type_commit(&parameters_struct);

	int root;
	MPI_Bcast(instance, 1, parameters_struct, root = 0, comm_shared);

	MPI_Type_free(&parameters_struct);
}

void initialize_problem(MPI_Comm comm_cart, instance_t * instance)
{
	const int split_dir = instance->local_subdomain_split_direction;

	const int LX = instance->subdomain_offsets[0] + instance->local_subdomain_offsets[0];
	const int LY = instance->subdomain_offsets[1] + instance->local_subdomain_offsets[1];
	const int LZ = instance->subdomain_offsets[2] + instance->local_subdomain_offsets[2];
	const int NX = instance->local_subdomain_sizes[0];
	const int NY = instance->local_subdomain_sizes[1];
	const int NZ = instance->local_subdomain_sizes[2];

	instance->F = (double*)calloc((size_t)NX * (size_t)NY * (size_t)NZ, sizeof(double));
	if (instance->F == NULL)
	{
		fprintf(stderr, "ERORR: memory allocation for F matrix failed! (Returned NULL)\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	double* F = instance->F;
	instance->dx[0] = 2.0 / (instance->domain_sizes[0] - 1.0);
	instance->dx[1] = 2.0 / (instance->domain_sizes[1] - 1.0);
	instance->dx[2] = 2.0 / (instance->domain_sizes[2] - 1.0);

	const double dx = instance->dx[0];
	const double dy = instance->dx[1];
	const double dz = instance->dx[2];
	const double m_alpha = -instance->alpha;

	double xval, yval, zval;

	for (int x = LX, i = 0; i < NX; x++, i++)
	{
		xval = -1.0 + dx * x;
		for (int y = LY, j = 0; j < NY; y++, j++)
		{
			yval = -1.0 + dy * y;
			for (int z = LZ, k = 0; k < NZ; z++, k++)
			{
				zval = -1.0 + dz * z;
				F[INDEX(i, j, k, NY, NZ)] =
					m_alpha * (1.0 - xval * xval) * (1.0 - yval * yval) * (1.0 - zval * zval)
					- 2.0 * (1.0 - yval * yval) * (1.0 - zval * zval)
					- 2.0 * (1.0 - xval * xval) * (1.0 - zval * zval)
					- 2.0 * (1.0 - xval * xval) * (1.0 - yval * yval);
			}
		}
	}
}

void close_problem(instance_t * instance)
{
	MPI_Win_fence(0, instance->win_U);
	MPI_Win_free(&instance->win_U);

	MPI_Win_fence(0, instance->win_Unew);
	MPI_Win_free(&instance->win_Unew);

	if (instance->F)
		free(instance->F);
}

void print_subdomain(double* mat, instance_t * instance, char* format)
{
	for (int i = 1; i <= instance->subdomain_sizes[0]; i++)
	{
		for (int j = 1; j <= instance->subdomain_sizes[1]; j++)
		{
			for (int k = 1; k <= instance->subdomain_sizes[2]; k++)
				printf(format,
					mat[INDEX(i, j, k,
						instance->subdomain_sizes[1] + 2,
						instance->subdomain_sizes[2] + 2)]);
			printf("\n");
		}
		printf("\n");
	}
}

void print_F(instance_t * instance, char* format)
{
	for (int i = 0; i < instance->local_subdomain_sizes[0]; i++)
	{
		for (int j = 0; j < instance->local_subdomain_sizes[1]; j++)
		{
			for (int k = 0; k < instance->local_subdomain_sizes[2]; k++)
				printf(format,
					instance->F[INDEX(i, j, k,
						instance->local_subdomain_sizes[1],
						instance->local_subdomain_sizes[2])]);
			printf("\n");
		}
		printf("\n");
	}
}

void setup_shared_and_heads(instance_t * instance, MPI_Comm * comm_shared, MPI_Comm * comm_head)
{
	int rank_world, rank_shared, nprocs_shared, color, root;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
	MPI_Bcast(&instance->use_shared_memory, 1, MPI_INT, root = 0, MPI_COMM_WORLD);
	MPI_Bcast(&instance->heads_per_shared_region, 1, MPI_INT, root = 0, MPI_COMM_WORLD);

	if (instance->use_shared_memory)
	{
		MPI_Comm comm_physical_shared;
		int rank_physical_shared, nprocs_physical_shared;
		MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_physical_shared);
		MPI_Comm_rank(comm_physical_shared, &rank_physical_shared);
		MPI_Comm_size(comm_physical_shared, &nprocs_physical_shared);

		/* 
			calculate a balanced split of the processes in the physical shared region.
			The reminder is distributed in a round-robin fashion among the firsts new islands.
			Example: nprocs_physical_shared = 8, heads_per_shared_region = 3.
			rank_physical_shared:	0,1,2,3,4,5,6,7
			color:					0,0,0,1,1,1,2,2
		*/
		if (instance->heads_per_shared_region > nprocs_physical_shared)
		{
			printf("WARNING: The requested number of heads per shared memory region (%i) exceeds "
				"the number of available processes in the shared region (%i).\n"
				"Setting 'Num heads per shared region' = %i.\n",
				instance->heads_per_shared_region,
				nprocs_physical_shared,
				nprocs_physical_shared);
			instance->heads_per_shared_region = nprocs_physical_shared;
			color = rank_physical_shared;
		}
		else
		{
			int stride, reminder;
			stride = nprocs_physical_shared / instance->heads_per_shared_region;
			reminder = nprocs_physical_shared % instance->heads_per_shared_region;
			if (rank_physical_shared / (stride + 1) < reminder)
				color = rank_physical_shared / (stride + 1);
			else
				color = (rank_physical_shared - reminder * (stride + 1)) / stride + reminder;
		}
		
		MPI_Comm_split(comm_physical_shared, color, 0, comm_shared);
	}
	else
	{
		color = rank_world;
		MPI_Comm_split(MPI_COMM_WORLD, color, 0, comm_shared);
	}

	MPI_Comm_rank(*comm_shared, &rank_shared);
	MPI_Comm_size(*comm_shared, &nprocs_shared);

	color = (rank_shared == 0) ? 1 : MPI_UNDEFINED;
	MPI_Comm_split(MPI_COMM_WORLD, color, 0, comm_head);
}

void setup_topology(MPI_Comm comm_head, int* nsplits_per_dim, int* coords, MPI_Comm * comm_cart)
{
	int rank_head, nprocs_head;
	if (comm_head != MPI_COMM_NULL)
	{
		MPI_Comm_rank(comm_head, &rank_head);
		MPI_Comm_size(comm_head, &nprocs_head);
	}
	*comm_cart = MPI_COMM_NULL;
	int rank_cart = MPI_PROC_NULL;
	int periods[DOMAIN_DIM];
	coords[0] = coords[1] = coords[2] = -1;
	memset(periods, 0x00, DOMAIN_DIM * sizeof(int));

	if (comm_head != MPI_COMM_NULL)
	{
		MPI_Dims_create(nprocs_head, DOMAIN_DIM, nsplits_per_dim);
		MPI_Cart_create(comm_head, DOMAIN_DIM, nsplits_per_dim, periods, 0, comm_cart);
	}
	if (*comm_cart != MPI_COMM_NULL)
	{
		MPI_Comm_rank(*comm_cart, &rank_cart);
		MPI_Cart_coords(*comm_cart, rank_cart, DOMAIN_DIM, coords);
	}
}

void compute_subdomains(MPI_Comm comm_head, int* coords, int* nsplits_per_dim, instance_t * instance)
{
	if (comm_head != MPI_COMM_NULL)
	{
		int divisor, reminder;
		for (int i = 0, index; i < DOMAIN_DIM; i++)
		{
			divisor = instance->domain_sizes[i] / nsplits_per_dim[i];
			reminder = instance->domain_sizes[i] % nsplits_per_dim[i];
			index = coords[i];
			instance->subdomain_offsets[i] = index * divisor + ((index <= reminder) ? index : reminder);
			index = coords[i] + 1;
			instance->subdomain_sizes[i] =
				index * divisor + ((index <= reminder) ? index : reminder) -
				instance->subdomain_offsets[i];
		}
	}
}

void compute_local_workload(MPI_Comm comm_shared, instance_t * instance)
{
	int rank_shared, nprocs_shared;
	MPI_Comm_rank(comm_shared, &rank_shared);
	MPI_Comm_size(comm_shared, &nprocs_shared);

	const int split_direction = instance->local_subdomain_split_direction;

	memset(instance->local_subdomain_offsets, 0x00, DOMAIN_DIM * sizeof(int));
	memcpy(instance->local_subdomain_sizes, instance->subdomain_sizes, DOMAIN_DIM * sizeof(int));

	int divisor, reminder, index;
	divisor = instance->subdomain_sizes[split_direction] / nprocs_shared;
	reminder = instance->subdomain_sizes[split_direction] % nprocs_shared;

	index = rank_shared;
	instance->local_subdomain_offsets[split_direction] =
		index * divisor + ((index <= reminder) ? index : reminder);

	index = rank_shared + 1;
	instance->local_subdomain_sizes[split_direction] =
		index * divisor + ((index <= reminder) ? index : reminder) -
		instance->local_subdomain_offsets[split_direction];
}

void allocate_shared_resources(MPI_Comm comm_cart, MPI_Comm comm_shared, instance_t * instance)
{
	MPI_Aint shared_size = (MPI_Aint)
		((uint64_t)instance->subdomain_sizes[0] + 2LL) *
		((uint64_t)instance->subdomain_sizes[1] + 2LL) *
		((uint64_t)instance->subdomain_sizes[2] + 2LL) * sizeof(double);

	if (comm_cart == MPI_COMM_NULL)
		shared_size = (MPI_Aint)0;

	MPI_Win_allocate_shared(shared_size, sizeof(double), MPI_INFO_NULL, comm_shared, &instance->U, &instance->win_U);
	MPI_Win_allocate_shared(shared_size, sizeof(double), MPI_INFO_NULL, comm_shared, &instance->Unew, &instance->win_Unew);

	int head = 0;
	int disp_unit;
	MPI_Aint actual_size;

	MPI_Win_shared_query(instance->win_U, head, &actual_size, &disp_unit, &instance->U);
	MPI_Win_shared_query(instance->win_Unew, head, &actual_size, &disp_unit, &instance->Unew);
}