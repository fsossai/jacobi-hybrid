#ifndef KERNEL_H
#define KERNEL_H

#include <mpi.h>
#include "data.h"

void compute_jacobi(MPI_Comm comm_cart, instance_t* instance);

#endif // !KERNEL_H


