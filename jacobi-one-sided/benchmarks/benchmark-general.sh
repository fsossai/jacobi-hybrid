#!/bin/bash

if (( $# < 3 )); then
	echo "Usage: $0 <cluster_name> <mode> <scaling>"
	echo "mode:	intranode | internode | topology"
	exit 1
fi

cluster_name=$1
mode=$2
scaling=$3
echo "Nodes = $SLURM_JOB_NODELIST"
echo "Maximum number of MPI processes: $SLURM_NTASKS"
time_start=$(date +%s)

if [[ $mode =~ ^(intranode|internode) ]]; then
	if [[ $scaling =~ ^(strong|weak) ]]; then
		./benchmark-scaling.sh $cluster_name $SLURM_NTASKS_PER_NODE $SLURM_JOB_NUM_NODES $mode $scaling
	else
		echo "ERROR: <scaling> must be 'strong' or 'weak'"
		exit 1
	fi
elif [ "$mode" = "topology" ]; then
	./benchmark-topology.sh $cluster_name $SLURM_NTASKS
else
	echo "ERROR: <mode> must be 'intranode', 'internode' or 'topology'"
	exit 1
fi

time_stop=$(date +%s)
echo [*] Job execution time: $((time_stop - time_start)) seconds