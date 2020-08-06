#include "kernel.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void compute_jacobi(MPI_Comm comm_cart, MPI_Comm comm_shared, instance_t* instance)
{
	const int OFFSETS[DOMAIN_DIM] =
	{
		instance->local_subdomain_offsets[0] + 1,
		instance->local_subdomain_offsets[1] + 1,
		instance->local_subdomain_offsets[2] + 1
	};
	const int U_NX = instance->subdomain_sizes[0] + 2;
	const int U_NY = instance->subdomain_sizes[1] + 2;
	const int U_NZ = instance->subdomain_sizes[2] + 2;

	double* U = instance->U;
	double* Unew = instance->Unew;

	const double* F = instance->F;
	const double ax = 1.0 / (instance->dx[0] * instance->dx[0]);
	const double ay = 1.0 / (instance->dx[1] * instance->dx[1]);
	const double az = 1.0 / (instance->dx[2] * instance->dx[2]);
	const double bb = -2.0 * (ax + ay + az) - instance->alpha;
	const double relax = instance->relaxation;

	MPI_Datatype facets_low_data[DOMAIN_DIM];
	MPI_Datatype facets_low_halo[DOMAIN_DIM];
	MPI_Datatype facets_high_data[DOMAIN_DIM];
	MPI_Datatype facets_high_halo[DOMAIN_DIM];
	int starts_low_data[DOMAIN_DIM];
	int starts_low_halo[DOMAIN_DIM];
	int starts_high_data[DOMAIN_DIM];
	int starts_high_halo[DOMAIN_DIM];
	int sizes[DOMAIN_DIM], subsizes[DOMAIN_DIM];
	memcpy(subsizes, instance->subdomain_sizes, DOMAIN_DIM * sizeof(int));

	for (int i = 0; i < DOMAIN_DIM; i++)
	{
		sizes[i] = instance->subdomain_sizes[i] + 2;
		starts_low_data[i] = starts_low_halo[i] = starts_high_data[i] = starts_high_halo[i] = 1;
	}
	for (int i = 0; i < DOMAIN_DIM; i++)
	{
		int temp = subsizes[i];
		subsizes[i] = 1;
		starts_low_halo[i] = 0;
		starts_high_data[i] = instance->subdomain_sizes[i];
		starts_high_halo[i] = instance->subdomain_sizes[i] + 1;

		MPI_Type_create_subarray(DOMAIN_DIM, sizes, subsizes, starts_low_data, MPI_ORDER_C, MPI_DOUBLE, &facets_low_data[i]);
		MPI_Type_create_subarray(DOMAIN_DIM, sizes, subsizes, starts_low_halo, MPI_ORDER_C, MPI_DOUBLE, &facets_low_halo[i]);
		MPI_Type_create_subarray(DOMAIN_DIM, sizes, subsizes, starts_high_data, MPI_ORDER_C, MPI_DOUBLE, &facets_high_data[i]);
		MPI_Type_create_subarray(DOMAIN_DIM, sizes, subsizes, starts_high_halo, MPI_ORDER_C, MPI_DOUBLE, &facets_high_halo[i]);
		MPI_Type_commit(&facets_low_data[i]);
		MPI_Type_commit(&facets_low_halo[i]);
		MPI_Type_commit(&facets_high_data[i]);
		MPI_Type_commit(&facets_high_halo[i]);

		// restoring initial general values
		subsizes[i] = temp;
		starts_low_halo[i] = starts_high_data[i] = starts_high_halo[i] = 1;
	}

	int rank_source[DOMAIN_DIM], rank_dest[DOMAIN_DIM];
	if (comm_cart != MPI_COMM_NULL)
	{
		for (int i = 0; i < DOMAIN_DIM; i++)
			MPI_Cart_shift(comm_cart, i, 1, &rank_source[i], &rank_dest[i]);
	}

	const int split_dir = instance->local_subdomain_split_direction;
	const int N[DOMAIN_DIM] =
	{ 
		instance->local_subdomain_sizes[0],
		instance->local_subdomain_sizes[1],
		instance->local_subdomain_sizes[2]
	};

	double residual, partial;
	MPI_Request requests[DOMAIN_DIM * 4];
	instance->performed_iterations = 0;
	double timer = -MPI_Wtime();
	for (int iteration = 0; iteration < instance->max_iterations; iteration++)
	{
		// halo exchange
		if (comm_cart != MPI_COMM_NULL)
		{
			int nreq = 0;
			for (int i = 0; i < DOMAIN_DIM; i++)
			{
				if (rank_source[i] != MPI_PROC_NULL)
				{
					MPI_Irecv(U, 1, facets_low_halo[i], rank_source[i], 0, comm_cart, &requests[nreq++]);
					MPI_Isend(U, 1, facets_low_data[i], rank_source[i], 0, comm_cart, &requests[nreq++]);
				}
				if (rank_dest[i] != MPI_PROC_NULL)
				{
					MPI_Isend(U, 1, facets_high_data[i], rank_dest[i], 0, comm_cart, &requests[nreq++]);
					MPI_Irecv(U, 1, facets_high_halo[i], rank_dest[i], 0, comm_cart, &requests[nreq++]);
				}
			}
			MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);
		}
		MPI_Barrier(comm_shared);

		// computation
		residual = 0.0;
		
		for (int i = 0, U_i = OFFSETS[0]; i < N[0]; i++, U_i++)
		{
			for (int j = 0, U_j = OFFSETS[1]; j < N[1]; j++, U_j++)
			{
				for (int k = 0, U_k = OFFSETS[2]; k < N[2]; k++, U_k++)
				{
					partial = (
						ax * (U[INDEX(U_i - 1, U_j, U_k, U_NY, U_NZ)] +
							U[INDEX(U_i + 1, U_j, U_k, U_NY, U_NZ)]) +
						ay * (U[INDEX(U_i, U_j - 1, U_k, U_NY, U_NZ)] +
							U[INDEX(U_i, U_j + 1, U_k, U_NY, U_NZ)]) +
						az * (U[INDEX(U_i, U_j, U_k - 1, U_NY, U_NZ)] +
							U[INDEX(U_i, U_j, U_k + 1, U_NY, U_NZ)]) +
						bb * U[INDEX(U_i, U_j, U_k, U_NY, U_NZ)] -
						F[INDEX(i, j, k, N[1], N[2])]
						) / bb;

					Unew[INDEX(U_i, U_j, U_k, U_NY, U_NZ)] =
						U[INDEX(U_i, U_j, U_k, U_NY, U_NZ)] - relax * partial;

					residual += partial * partial;
				}
			}
		}
		// getting total residual inside a shared memory island
		double local_residual;
		MPI_Reduce(&residual, &local_residual, 1, MPI_DOUBLE, MPI_SUM, 0, comm_shared);

		// getting total residual of the whole computation
		double total_residual;
		if (comm_cart != MPI_COMM_NULL)
			MPI_Allreduce(&local_residual, &total_residual, 1, MPI_DOUBLE, MPI_SUM, comm_cart);

		total_residual = sqrt(total_residual) /
			instance->subdomain_sizes[0] *
			instance->subdomain_sizes[1] *
			instance->subdomain_sizes[2];

		// swapping pointers
		double* temp = U;
		instance->U = U = Unew;
		Unew = temp;

		instance->residual = total_residual;
		instance->performed_iterations++;
	}
	timer += MPI_Wtime();
	instance->total_computation_time = timer;

	for (int i = 0; i < DOMAIN_DIM; i++)
	{
		MPI_Type_free(&facets_low_data[i]);
		MPI_Type_free(&facets_low_halo[i]);
		MPI_Type_free(&facets_high_data[i]);
		MPI_Type_free(&facets_high_halo[i]);
	}
}
