#include <mpi.h>
#include <stddef.h>
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