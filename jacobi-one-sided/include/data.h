#ifndef DATA_H
#define DATA_H

#include <stdio.h>

#define DOMAIN_DIM 3
#define INDEX3D(i,j,k,N2,N3) ((i)*(N2)*(N3) + (j)*(N3) + (k))
#define INDEX2D(i,j,N2) ((i)*(N2) (j))

typedef struct
{
	int domain_sizes[DOMAIN_DIM];
	int subdomain_sizes[DOMAIN_DIM];
	int subdomain_offsets[DOMAIN_DIM];

	double alpha;
	double residual;
	double relaxation;
	double tolerance;
	int max_iterations;

	double* U;
	double* F;

	double dx[DOMAIN_DIM];

} instance_t;

void read_input(FILE* stream, instance_t* instance);
void broadcast_input_data(instance_t* instance);
void initialize_problem(instance_t* instance);
void close_problem(instance_t* instance);
void print_subdomain(double* mat, instance_t* instance, char* format);

#endif // !DATA_H