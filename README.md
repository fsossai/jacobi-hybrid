# Hybrid implementation of Jacobi algorithm
(Work in progress!)

## Summary
This is an MPI parallel implementation of the Jacobi algorithm that is meant to exploit both distributed and shared memory.
This does not make use of threads: physical shared memory is managed through MPI one-sided communication features of MPI 3.0.

## Motivation and Description
This project aims to develop a numerical solver of a differential equation using the so called Jacobi algorithm. In layman's terms this simple algorithm iteratively computes a grid of points from a previous one in which every point depends only on its neighbours in the previous grid. If the grid gets splitted in sub-grids (in order to feed different distributed nodes), points at the boundaries will depend on points which are not in the same sub-grid. A data communication among sub-grids is therefore needed; this aspect will mostly characterise the performance of the program.
Jacobi is not the best candidate to carry out the numerical job, but here the main goal is to deal with all the issues that come with the distributed + shared memory paradigm and providing an implementation suitable for supercomputers.

## Tools
Now let's talk about nuts and bolts. Some hybrid programming models combine two different standards like MPI and OpenMP to take advantage of the underlying architecture, however in my project I will pursue different path, using MPI 3.0 shared memory that allows to exploit shared and distributed memory at the same time by means of the so called one-sided communications feature. 
