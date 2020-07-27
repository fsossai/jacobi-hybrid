#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>
#include "data.h"

void read_input(FILE* stream, instance_t *instance)
{
	fscanf(stream, "%i", &instance->domain_x_size);
	fscanf(stream, "%i", &instance->domain_y_size);
	fscanf(stream, "%i", &instance->domain_z_size);
	fscanf(stream, "%lf", &instance->alpha);
	fscanf(stream, "%lf", &instance->relaxation);
	fscanf(stream, "%lf", &instance->tolerance);
	fscanf(stream, "%i", &instance->max_iterations);

	#ifdef _DEBUG
	printf("dims\t\t: (%i,%i,%i)\n", instance->domain_x_size, instance->domain_y_size, instance->domain_z_size);
	printf("alpha\t\t: %lf\n", instance->alpha);
	printf("relaxation\t: %lf\n", instance->relaxation);
	printf("tolerance\t: %lf\n", instance->tolerance);
	printf("iterations\t: %i\n", instance->max_iterations);
	#endif // DEBUG

}

void broadcast_input_data(instance_t *instance)
{
	MPI_Aint displacements[7];
	displacements[0] = (MPI_Aint)offsetof(instance_t, domain_x_size);
	displacements[1] = (MPI_Aint)offsetof(instance_t, domain_y_size);
	displacements[2] = (MPI_Aint)offsetof(instance_t, domain_z_size);
	displacements[3] = (MPI_Aint)offsetof(instance_t, alpha);
	displacements[4] = (MPI_Aint)offsetof(instance_t, relaxation);
	displacements[5] = (MPI_Aint)offsetof(instance_t, tolerance);
	displacements[6] = (MPI_Aint)offsetof(instance_t, max_iterations);

	int block_lengths[7];
	for (int i = 0; i < 7; i++)
		block_lengths[i] = 1;

	MPI_Datatype types[7] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };

	MPI_Datatype parameters_struct;
	MPI_Type_create_struct(7, block_lengths, displacements, types, &parameters_struct);
	MPI_Type_commit(&parameters_struct);

	MPI_Bcast(instance, 1, parameters_struct, 0, MPI_COMM_WORLD);

	MPI_Type_free(&parameters_struct);
}

void initialize_problem(instance_t *instance)
{
	const int LX = instance->offset_x;
	const int LY = instance->offset_y;
	const int LZ = instance->offset_z;
	const int NX = instance->subdomain_x_size;
	const int NY = instance->subdomain_y_size;
	const int NZ = instance->subdomain_z_size;
	const int N = NX * NY * NZ;
	
	instance->U = (double*)calloc(N, sizeof(double));
	instance->F = (double*)malloc(N * sizeof(double));
	instance->dx = 2.0 / (instance->domain_x_size - 1);
	instance->dy = 2.0 / (instance->domain_y_size - 1);
	instance->dz = 2.0 / (instance->domain_z_size - 1);

	double* U = instance->U;
	double* F = instance->F;
	const double dx = instance->dx;
	const double dy = instance->dy;
	const double dz = instance->dz;
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
				F[INDEX(i, j, k, NY, NZ)] =
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
	for (int i = 0; i<instance->subdomain_x_size; i++)
	{
		for (int j = 0; j < instance->subdomain_y_size; j++)
		{
			for (int k = 0; k < instance->subdomain_z_size; k++)
				printf(format,
					mat[INDEX(i, j, k,
						instance->subdomain_y_size,
						instance->subdomain_z_size)]);
			printf("\n");
		}
		printf("\n");
	}
}