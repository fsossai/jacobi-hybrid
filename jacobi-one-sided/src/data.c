#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>
#include "data.h"

void read_input(FILE* stream, instance_t *instance)
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

void broadcast_input_data(instance_t *instance)
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

	MPI_Bcast(instance, 1, parameters_struct, 0, MPI_COMM_WORLD);

	MPI_Type_free(&parameters_struct);
}

void initialize_problem(instance_t *instance)
{
	const int LX = instance->subdomain_offsets[0];
	const int LY = instance->subdomain_offsets[1];
	const int LZ = instance->subdomain_offsets[2];
	const int NX = instance->subdomain_sizes[0];
	const int NY = instance->subdomain_sizes[1];
	const int NZ = instance->subdomain_sizes[2];
	const int N = NX * NY * NZ;
	
	instance->U = (double*)calloc(N, sizeof(double));
	instance->F = (double*)malloc(N * sizeof(double));
	instance->dx[0] = 2.0 / (instance->domain_sizes[0] - 1);
	instance->dx[1] = 2.0 / (instance->domain_sizes[1] - 1);
	instance->dx[2] = 2.0 / (instance->domain_sizes[2] - 1);

	double* U = instance->U;
	double* F = instance->F;
	const double dx = instance->dx[0];
	const double dy = instance->dx[1];
	const double dz = instance->dx[2];
	const double m_alpha = - instance->alpha;

	double xval, yval, zval;

	for (int x = LX, i=0; i<NX; x++, i++)
	{
		xval = -1.0 + dx * x;
		for (int y = LY, j=0; j<NY; y++, j++)
		{
			yval = -1.0 + dy * y;
			for (int z = LZ, k=0; k<NZ; z++, k++)
			{
				zval = -1.0 + dz * z;
				F[INDEX3D(i, j, k, NY, NZ)] =
					m_alpha * (1.0 - xval * xval) * (1.0 - yval * yval) * (1.0 - zval * zval) +
					2.0 * (-2.0 + xval * xval + yval * yval + zval * zval);
			}
		}
	}
}

void close_problem(instance_t* instance)
{
	if (instance->U)
		free(instance->U);
	if (instance->F)
		free(instance->F);
}

void print_subdomain(double* mat, instance_t *instance, char *format)
{
	for (int i = 0; i<instance->subdomain_sizes[0]; i++)
	{
		for (int j = 0; j < instance->subdomain_sizes[1]; j++)
		{
			for (int k = 0; k < instance->subdomain_sizes[2]; k++)
				printf(format,
					mat[INDEX3D(i, j, k,
						instance->subdomain_sizes[1],
						instance->subdomain_sizes[2])]);
			printf("\n");
		}
		printf("\n");
	}
}