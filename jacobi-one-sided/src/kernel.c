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
	double* F = instance->F;
	double* Unew = (double*)malloc(U_NX * U_NY * U_NZ * sizeof(double));
	const double dx = instance->dx[0];
	const double dy = instance->dx[1];
	const double dz = instance->dx[2];
	const double m_alpha = -instance->alpha;

	MPI_Datatype facets[DOMAIN_DIM];
	int starts[DOMAIN_DIM], sizes[DOMAIN_DIM], subsizes[DOMAIN_DIM];
	memset(starts, 0x00, DOMAIN_DIM * sizeof(int));
	memcpy(subsizes, instance->subdomain_sizes, DOMAIN_DIM * sizeof(int));

	for (int i = 0; i < DOMAIN_DIM; i++)
		sizes[i] = instance->subdomain_sizes[i] + 2;
	for (int i = 0; i < DOMAIN_DIM; i++)
	{
		int temp = subsizes[i];
		subsizes[i] = 1;
		MPI_Type_create_subarray(DOMAIN_DIM, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &facets[i]);
		MPI_Type_commit(&facets[i]);
		subsizes[i] = temp;
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
			printf("i:%i source:%i, dest:%i\n", i, rank_source[i], rank_dest[i]); fflush(stdout); //D
			if (rank_source[i] != MPI_PROC_NULL)
			{
				printf("[%i] Low comm\n",i);
				MPI_Irecv(&U[INDEX3D(i != 0, i != 1, i != 2, U_NY, U_NZ)], 1, facets[i], rank_source[i], 0, comm_cart, &requests[nreq++]);
				MPI_Isend(&U[INDEX3D(1, 1, 1, U_NY, U_NZ)], 1, facets[i], rank_source[i], 0, comm_cart, &requests[nreq++]);
			}
			if (rank_dest[i] != MPI_PROC_NULL)
			{
				printf("[%i] High comm\n", i);
				MPI_Isend(&U[INDEX3D(
					(i == 0) ? N[i] : 1,
					(i == 1) ? N[i] : 1,
					(i == 2) ? N[i] : 1,
					U_NY, U_NZ)], 1, facets[i], rank_dest[i], 0, comm_cart, &requests[nreq++]);
				MPI_Irecv(&U[INDEX3D(
					(i == 0) ? (1 + N[i]) : 1,
					(i == 1) ? (1 + N[i]) : 1,
					(i == 2) ? (1 + N[i]) : 1,
					U_NY, U_NZ)], 1, facets[i], rank_dest[i], 0, comm_cart, &requests[nreq++]);
			}
		}
		MPI_Waitall(nreq, requests, MPI_STATUSES_IGNORE);

		// computation
		for (int i = 1; i<=N[0]; i++)
		{
			for (int j = 1; j<=N[1]; j++)
			{
				for (int k = 1; k<=N[2]; k++)
				{
					Unew[INDEX3D(i, j, k, U_NY, U_NZ)] = (
						U[INDEX3D(i - 1, j, k, U_NY, U_NZ)] +
						U[INDEX3D(i, j - 1, k, U_NY, U_NZ)] +
						U[INDEX3D(i, j, k - 1, U_NY, U_NZ)] +
						U[INDEX3D(i, j, k + 1, U_NY, U_NZ)] +
						U[INDEX3D(i, j + 1, k, U_NY, U_NZ)] +
						U[INDEX3D(i + 1, j, k, U_NY, U_NZ)]
						) / 6.0;
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
