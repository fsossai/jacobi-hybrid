#ifndef KERNEL_H
#define KERNEL_H

#include <mpi.h>
#include "data.h"

#define KERNEL_FLOPS 13
#define RESIDUAL_CHECK 100

void compute_jacobi(MPI_Comm comm_cart, MPI_Comm comm_shared, instance_t* instance);

#endif // !KERNEL_H


