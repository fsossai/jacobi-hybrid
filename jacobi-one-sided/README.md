## Compiling

In this repository, the two version of the Jacobi have been splitted in two different branches.
The _master_ branch contains the 3D version and the benchmarking infrastructure, whereas
_2d-version_ contains only the 2D version.

To compile everything:
Remember to load the necessary modules if not already loaded, e.g.:
```
module load intel/18 intel-mpi/2018
```
Suppose git is set to the master branch. Then compile first the 3D version and then 2D:
```
make release
git checkout 2d-version
make release
```
`release` will use `-xHost`, generating instructions for the highest instruction set available on the
compilation host processor. Use `make skylake` or `make ivybridge` instead, if you what to specifically
generate executables for one of those two microarchitectures.

For **debugging** purposes, it is possible to show all informations about the domain decompositions,
cartesian topology, and communicators in general during the execution of the program.
Use `-DDEBUG` as an additional argument to the compiler command line, or simply compile with `make debug`.
