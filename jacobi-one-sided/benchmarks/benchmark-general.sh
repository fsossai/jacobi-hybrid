#!/bin/bash

if (( $# < 2 )); then
	echo "Usage: $0 <cluster_name> <mode>"
	echo "mode:	intranode | internode | topology"
	exit 1
fi

cluster_name=$1
mode=$2
echo "Nodes = $SLURM_JOB_NODELIST"
echo "Maximum number of MPI processes: $SLURM_NTASKS"
time_start=$(date +%s)

if [[ $mode =~ ^(intranode|internode) ]]; then
	./benchmark-scaling.sh $cluster_name $SLURM_NTASKS_PER_NODE $SLURM_JOB_NUM_NODES $mode
elif [ "$mode" = "topology" ]; then
	./benchmark-topology.sh $cluster_name $SLURM_NTASKS
else
	echo "ERROR: Unrecognized argument '$mode'"
	echo "Available arguments: 'intranode', 'internode', 'topology'"
	exit 1
fi

time_stop=$(date +%s)
echo [*] Job execution time: $((time_stop - time_start)) seconds