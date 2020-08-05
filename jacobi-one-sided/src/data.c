#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"

void read_input(FILE* stream, instance_t* instance)
{
	for (int i = 0; i < DOMAIN_DIM; i++)
		fscanf(stream, "%i", &instance->domain_sizes[i]);
	fscanf(stream, "%lf", &instance->alpha);
	fscanf(stream, "%lf", &instance->relaxation);
	fscanf(stream, "%lf", &instance->tolerance);
	fscanf(stream, "%i", &instance->max_iterations);

	#ifdef _DEBUG
	printf("dims\t\t: (%i,%i,%i)\n",
		instance->domain_sizes[0], instance->domain_sizes[1], instance->domain_sizes[2]);
	printf("alpha\t\t: %lf\n", instance->alpha);
	printf("relaxation\t: %lf\n", instance->relaxation);
	printf("tolerance\t: %lf\n", instance->tolerance);
	printf("iterations\t: %i\n", instance->max_iterations);
	#endif // _DEBUG

}

void broadcast_input_data_head(MPI_Comm comm_head, instance_t * instance)
{
	MPI_Aint displacements[5];
	displacements[0] = (MPI_Aint)offsetof(instance_t, domain_sizes);
	displacements[1] = (MPI_Aint)offsetof(instance_t, alpha);
	displacements[2] = (MPI_Aint)offsetof(instance_t, relaxation);
	displacements[3] = (MPI_Aint)offsetof(instance_t, tolerance);
	displacements[4] = (MPI_Aint)offsetof(instance_t, max_iterations);

	int block_lengths[5];
	block_lengths[0] = DOMAIN_DIM;
	for (int i = 1; i < 5; i++)
		block_lengths[i] = 1;

	MPI_Datatype types[5] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };

	MPI_Datatype parameters_struct;
	MPI_Type_create_struct(5, block_lengths, displacements, types, &parameters_struct);
	MPI_Type_commit(&parameters_struct);

	int root;
	if (comm_head != MPI_COMM_NULL)
		MPI_Bcast(instance, 1, parameters_struct, root = 0, comm_head);

	MPI_Type_free(&parameters_struct);
}

void broadcast_data_shared(MPI_Comm comm_shared, instance_t * instance)
{
	MPI_Aint displacements[5];
	displacements[0] = (MPI_Aint)offsetof(instance_t, subdomain_sizes);
	displacements[1] = (MPI_Aint)offsetof(instance_t, alpha);
	displacements[2] = (MPI_Aint)offsetof(instance_t, relaxation);
	displacements[3] = (MPI_Aint)offsetof(instance_t, tolerance);
	displacements[4] = (MPI_Aint)offsetof(instance_t, max_iterations);

	int block_lengths[5];
	block_lengths[0] = DOMAIN_DIM;
	for (int i = 1; i < 5; i++)
		block_lengths[i] = 1;

	MPI_Datatype types[5] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };

	MPI_Datatype parameters_struct;
	MPI_Type_create_struct(5, block_lengths, displacements, types, &parameters_struct);
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

	instance->U = NULL;
	instance->F = (double*)malloc(NX * NY * NZ * sizeof(double));

	if (comm_cart != MPI_COMM_NULL)
	{
		instance->U = (double*)calloc(
			(instance->subdomain_sizes[0] + 2) *
			(instance->subdomain_sizes[1] + 2) *
			(instance->subdomain_sizes[2] + 2), sizeof(double));
	}
	
	double* F = instance->F;
	instance->dx[0] = 2.0 / (instance->domain_sizes[0] - 1);
	instance->dx[1] = 2.0 / (instance->domain_sizes[1] - 1);
	instance->dx[2] = 2.0 / (instance->domain_sizes[2] - 1);

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
	if (instance->U)
		free(instance->U);
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

void setup_shared_and_heads(int nheads_per_node, MPI_Comm * comm_shared, MPI_Comm * comm_head)
{
	int rank_shared, nprocs_shared, color;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_shared);
	MPI_Comm_rank(*comm_shared, &rank_shared);
	MPI_Comm_size(*comm_shared, &nprocs_shared);

	int stride = nprocs_shared / nheads_per_node;
	if (stride > 0)
		color = (rank_shared % stride == 0 &&
			rank_shared / stride < nheads_per_node) ? 1 : MPI_UNDEFINED;
	else
		color = 1;

	MPI_Comm_split(MPI_COMM_WORLD, color, 0, comm_head);
}

void setup_topology(MPI_Comm comm_head, int* nsplits_per_dim, int* coords, MPI_Comm * comm_cart)
{
	int rank_head, nprocs_head = 0;
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
	memset(nsplits_per_dim, 0x00, DOMAIN_DIM * sizeof(int));

	#ifdef NO_SPLIT_X
	nsplits_per_dim[0] = 1;
	#endif
	#ifdef NO_SPLIT_Y
	nsplits_per_dim[1] = 1;
	#endif
	#ifdef NO_SPLIT_Z
	nsplits_per_dim[2] = 1;
	#endif

	if (nprocs_head > 0)
		MPI_Dims_create(nprocs_head, DOMAIN_DIM, nsplits_per_dim);
	if (comm_head != MPI_COMM_NULL)
		MPI_Cart_create(comm_head, DOMAIN_DIM, nsplits_per_dim, periods, 0, comm_cart);
	if (*comm_cart != MPI_COMM_NULL)
	{
		MPI_Comm_rank(*comm_cart, &rank_cart);
		MPI_Cart_coords(*comm_cart, rank_cart, DOMAIN_DIM, coords);
	}
}

void compute_subdomains(MPI_Comm comm_cart, int* coords, int* nsplits_per_dim, instance_t * instance)
{
	if (comm_cart != MPI_COMM_NULL)
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

	int split_direction = DEFAULT_LOCAL_SPLIT;
	#ifdef LOCAL_SPLIT_X
	split_direction = 0;
	#endif
	#ifdef LOCAL_SPLIT_Y
	split_direction = 1;
	#endif
	#ifdef LOCAL_SPLIT_Z
	split_direction = 2;
	#endif
	instance->local_subdomain_split_direction = split_direction;

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