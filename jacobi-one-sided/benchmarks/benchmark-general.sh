#!/bin/bash

if (( $# < 1 )); then
	echo "Usage: $0 <mode>"
	echo "mode:	intranode | internode | topology"
	exit 1
fi

cluster_name=vsc3		# Just to distinguish the slurm output files easily
hpn=2					# Number of heads per node (valid only for 'topology' option)
echo "Nodes = $SLURM_JOB_NODELIST"
echo "Maximum number of MPI processes: $SLURM_NTASKS"
mode=$1
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