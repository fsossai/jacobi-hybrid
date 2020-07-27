#ifndef DATA_H
#define DATA_H

#include <stdio.h>

#define INDEX(i,j,k,N1,N2) ((i)*(N1)*(N2) + (j)*(N2) + (k))

typedef struct
{
	int domain_x_size;
	int domain_y_size;
	int domain_z_size;

	int subdomain_x_size;
	int subdomain_y_size;
	int subdomain_z_size;

	int offset_x;
	int offset_y;
	int offset_z;

	double alpha;
	double residual;
	double relaxation;
	double tolerance;
	int max_iterations;

	double* U;
	double* F;

	double dx, dy, dz;

} instance_t;

void read_input(FILE* stream, instance_t* instance);
void broadcast_input_data(instance_t* instance);
void initialize_problem(instance_t* instance);
void close_problem(instance_t* instance);
void print_subdomain(double* mat, instance_t* instance, char* format);

#endif // !DATA_H