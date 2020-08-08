#ifndef DATA_H
#define DATA_H

#include <stdio.h>
#include <mpi.h>

#define DOMAIN_DIM 3
#define INDEX(i,j,k,N2,N3) ((i)*(N2)*(N3) + (j)*(N3) + (k))
#define DEFAULT_SHARED_SPLIT_DIRECTION 0
#define DEFAULT_USE_SHARED_MEMORY 1
#define DEFAULT_HEADS_PER_SHARED_REGION 1
#define STRING_FILE_NAME_MAX_SIZE 128

typedef struct
{
	int domain_sizes[DOMAIN_DIM];
	int subdomain_sizes[DOMAIN_DIM];
	int subdomain_offsets[DOMAIN_DIM];
	int local_subdomain_split_direction;
	int local_subdomain_offsets[DOMAIN_DIM];
	int local_subdomain_sizes[DOMAIN_DIM];
	int cart_splits[DOMAIN_DIM];
	int use_shared_memory;
	int heads_per_shared_region;

	char input_file_name[STRING_FILE_NAME_MAX_SIZE];
	FILE* input_stream;
	double alpha;
	double residual;
	double relaxation;
	double tolerance;
	double dx[DOMAIN_DIM];
	int max_iterations;
	int performed_iterations;
	double total_computation_time;

	double* U;
	double* Unew;
	double* F;

	MPI_Win win_U;
	MPI_Win win_Unew;

} instance_t;

void read_input(instance_t* instance);
void broadcast_input_data_head(MPI_Comm comm_head, instance_t* instance);
void broadcast_data_shared(MPI_Comm comm_shared, instance_t* instance);
void initialize_problem(MPI_Comm comm_cart, instance_t* instance);
void close_problem(instance_t* instance);
void print_subdomain(double* mat, instance_t* instance, char* format);
void print_F(instance_t* instance, char* format);
void setup_shared_and_heads(instance_t* instance, MPI_Comm* comm_shared, MPI_Comm* comm_head);
void setup_topology(MPI_Comm comm_head, int* nsplits_per_dim, int* coords, MPI_Comm* comm_cart);
void compute_subdomains(MPI_Comm comm_head, int* coords, int* nsplits_per_dim, instance_t* instance);
void compute_local_workload(MPI_Comm comm_shared, instance_t* instance);
void allocate_shared_resources(MPI_Comm comm_cart, MPI_Comm comm_shared, instance_t* instance);

#endif // !DATA_H