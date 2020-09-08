This folder contains only job files and scripts to extensively benchmark either the
pure and the hybrid version of the Jacobi.

## Overview
- **intranode**: Tests the scalability using, one by one, the processing units of a single node.
The baseline for parallel efficiency is a single MPI process.
- **internode**: Tests the scalability using, one by one, a given number of computing nodes in which all available cores
are used. The baseline for parallel efficiency is an entire node.
`benchmark-vsc3.job`, `benchmark-vsc3plus.job`, and `benchmark-vsc4.job` are job files already configured
for the three different architectures.
Both _intranode_ and _internode_ types of scaling can be benchmarked using `benchmark-general.sh` which is
the entry point for the two and for the _topology_ scaling as well.
By default, intranode and internode modes will produce results for both _strong_ and _weak_ scaling.

Parameters like the problem size and/or number of iterations to be used during the intranode/internode scaling benchmarks
can be changed in `benchmark-scaling.sh`. The latter will automatically test the 2D and 3D versions.

### Call structure

This shell scripts are thought to be modular and flexible.

```
|-> benchmark-general.sh
        |-> benchmark-scaling.sh
                |-> intranode-scaling.sh
                |-> internode-scaling.sh
        |-> benchmark-topology.sh
                |-> topology-scaling.sh 
```

### Examples
To submit to Slurm an internode scaling benchmark on VSC4:
```sbatch benchmark-vsc4.job internode```

The job will produce a set of files in a csv and table format:
```
bm_vsc4_intra_strong2d_2020-09-08_19-21.csv
bm_vsc4_intra_strong2d_2020-09-08_19-21.table
bm_vsc4_intra_strong3d_2020-09-08_19-27.csv
bm_vsc4_intra_strong3d_2020-09-08_19-27.table
...
```

Running, `cat bm_vsc4*table` it is possible to show the results collected in the table-like format:
```
2D Intranode strong scaling on vsc4 (1 runs per method)
Instance  Procs  |  P:Topology  P:Efficiency  P:Performance  P:Time  |  H:Topology  H:Efficiency  H:Performance  H:Time
(9000)^2  1      |  1x1         1.000         2.762          19.058  |  1x1         1.000         2.767          19.022
(9000)^2  2      |  2x1         0.975         5.392          9.764   |  1x1         0.474         2.626          20.048
(9000)^2  3      |  3x1         0.962         7.980          6.597   |  1x1         0.473         3.935          13.377
(9000)^2  4      |  2x2         0.947         10.473         5.027   |  1x1         0.470         5.212          10.101
...
```