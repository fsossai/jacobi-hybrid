#ifndef DATA_H
#define DATA_H

#include <stdio.h>

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

} instance_t;

void read_input(FILE* stream, instance_t* instance);
void broadcast_input_data(instance_t* instance);

#endif // !DATA_H