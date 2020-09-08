#!/bin/bash

if (( $# < 4 )); then
    echo "Usage: $0 <cluster_name> <ppn> <nnodes> <mode>"
    echo "mode: intranode | internode"
    exit 1
fi

cluster_name=$1
ppn=$2
nnodes=$3
mode=$4

if [ "$mode" = "intranode" ]; then
	nnodes=""
elif [ "$mode" = "internode" ]; then
	nnodes=$SLURM_JOB_NUM_NODES
else
    echo "ERROR: Wrong argument for <mode> ($mode)"
    echo "mode: intranode | internode"
    exit 1
fi

echo "Running ${mode}-scaling.sh (2D / strong) on $cluster_name"
./${mode}-scaling.sh 2 44000 50 $ppn $nnodes strong $cluster_name

echo "Running ${mode}-scaling.sh (3D / strong) on $cluster_name"
./${mode}-scaling.sh 3 1200 50 $ppn $nnodes strong $cluster_name

echo "Running ${mode}-scaling.sh (2D / weak) on $cluster_name"
./${mode}-scaling.sh 2 44000 50 $ppn $nnodes weak $cluster_name

echo "Running ${mode}-scaling.sh (3D / weak) on $cluster_name"
./${mode}-scaling.sh 3 1200 50 $ppn $nnodes weak $cluster_name