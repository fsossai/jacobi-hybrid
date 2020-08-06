#ifndef DATA_H
#define DATA_H

#include <stdio.h>

#define DOMAIN_DIM 3
#define INDEX(i,j,k,N2,N3) ((i)*(N2)*(N3) + (j)*(N3) + (k))
#define DEFAULT_LOCAL_SPLIT 0
#define _CONFIRM { int __rank; MPI_Comm_rank(MPI_COMM_WORLD,&__rank); printf("[%2i] OK\n",__rank); fflush(stdout); }

typedef struct
{
	int domain_sizes[DOMAIN_DIM];
	int subdomain_sizes[DOMAIN_DIM];
	int subdomain_offsets[DOMAIN_DIM];
	int local_subdomain_split_direction;
	int local_subdomain_offsets[DOMAIN_DIM];
	int local_subdomain_sizes[DOMAIN_DIM];

	double alpha;
	double residual;
	double relaxation;
	double tolerance;
	int max_iterations;
	int performed_iterations;
	double total_computation_time;

	double* U;
	double* Unew;
	double* F;

	MPI_Win win_U;
	MPI_Win win_Unew;

	double dx[DOMAIN_DIM];

} instance_t;

void read_input(FILE* stream, instance_t* instance);
void broadcast_input_data_head(MPI_Comm comm_head, instance_t* instance);
void broadcast_data_shared(MPI_Comm comm_shared, instance_t* instance);
void initialize_problem(MPI_Comm comm_cart, instance_t* instance);
void close_problem(instance_t* instance);
void print_subdomain(double* mat, instance_t* instance, char* format);
void print_F(instance_t* instance, char* format);
void setup_shared_and_heads(int nheads_per_node, MPI_Comm* comm_shared, MPI_Comm* comm_head);
void setup_topology(MPI_Comm comm_head, int* nsplits_per_dim, int* coords, MPI_Comm* comm_cart);
void compute_subdomains(MPI_Comm comm_cart, int* coords, int* nsplits_per_dim, instance_t* instance);
void compute_local_workload(MPI_Comm comm_shared, instance_t* instance);
void allocate_shared_resources(MPI_Comm comm_cart, MPI_Comm comm_shared, instance_t* instance);

#endif // !DATA_H