#include "kernel.h"
#include <stdlib.h>
#include <string.h>

void compute_jacobi(MPI_Comm comm_cart, instance_t* instance)
{
	if (comm_cart == MPI_COMM_NULL)
		return;

	int N[DOMAIN_DIM];
	N[0] = instance->subdomain_sizes[0];
	N[1] = instance->subdomain_sizes[1];
	N[2] = instance->subdomain_sizes[2];
	const int U_NX = N[0] + 2;
	const int U_NY = N[1] + 2;
	const int U_NZ = N[2] + 2;

	double* U = instance->U;
	double* Unew = (double*)calloc(U_NX * U_NY * U_NZ, sizeof(double));
	const double* F = instance->F;
	//double* Unew = (double*)malloc(U_NX * U_NY * U_NZ * sizeof(double));
	const double ax = 1.0 / (instance->dx[0] * instance->dx[0]);
	const double ay = 1.0 / (instance->dx[1] * instance->dx[1]);
	const double az = 1.0 / (instance->dx[2] * instance->dx[2]);
	const double bb = -2.0 * (ax + ay + az) - instance->alpha;

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
		sizes[i] = N[i] + 2;
		starts_low_data[i] = starts_low_halo[i] = starts_high_data[i] = starts_high_halo[i] = 1;
	}
	for (int i = 0; i < DOMAIN_DIM; i++)
	{
		int temp = subsizes[i];
		subsizes[i] = 1;
		starts_low_halo[i] = 0;
		starts_high_data[i] = N[i];
		starts_high_halo[i] = N[i] + 1;

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
	for (int i = 0; i < DOMAIN_DIM; i++)
		MPI_Cart_shift(comm_cart, i, 1, &rank_source[i], &rank_dest[i]);

	MPI_Request requests[DOMAIN_DIM * 4];
	for (int iteration = 0; iteration < instance->max_iterations; iteration++)
	{
		// halo exchange
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

		// computation
		for (int i = 1; i <= N[0]; i++)
		{
			for (int j = 1; j <= N[1]; j++)
			{
				for (int k = 1; k <= N[2]; k++)
				{
					Unew[INDEX3D(i, j, k, U_NY, U_NZ)] = (
						ax * (U[INDEX3D(i - 1, j, k, U_NY, U_NZ)] +
							U[INDEX3D(i + 1, j, k, U_NY, U_NZ)]) +
						ay * (U[INDEX3D(i, j - 1, k, U_NY, U_NZ)] +
							U[INDEX3D(i, j + 1, k, U_NY, U_NZ)]) +
						az * (U[INDEX3D(i, j, k - 1, U_NY, U_NZ)] +
							U[INDEX3D(i, j, k + 1, U_NY, U_NZ)]) +
						bb * U[INDEX3D(i, j, k, U_NY, U_NZ)] -
						F[INDEX3D(i, j, k, U_NY, U_NZ)]
						) / bb;

					/*fLRes = (ax * (UOLD(j, i - 1) + UOLD(j, i + 1))
						+ ay * (UOLD(j - 1, i) + UOLD(j + 1, i))
						+ b * UOLD(j, i) - F(j, i)) / b;

					U(j, i) = UOLD(j, i) - data->fRelax * fLRes;

					residual += fLRes * fLRes;*/
				}
			}
		}

		// swapping pointers
		double* temp = U;
		instance->U = U = Unew;
		Unew = temp;
	}
	free(Unew);
}
