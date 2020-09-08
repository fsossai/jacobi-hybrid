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
	nnodes=$SLURM_JOB_NUM_NODES
elif [ "$mode" = "internode" ]; then
	nnodes=""
else
    echo "ERROR: Wromg argument for <mode> ($mode)"
    echo "mode: intranode | internode"
    exit 1
fi

echo "Running ${mode}-scaling.sh (2D / strong) on $cluster_name"
./${mode}-scaling.sh 44000 50 2 $ppn $nnodes strong $cluster_name

echo "Running ${mode}-scaling.sh (3D / strong) on $cluster_name"
./${mode}-scaling.sh 1200 50 3 $ppn $nnodes strong $cluster_name

echo "Running ${mode}-scaling.sh (2D / weak) on $cluster_name"
./${mode}-scaling.sh 44000 50 2 $ppn $nnodes weak $cluster_name

echo "Running ${mode}-scaling.sh (3D / weak) on $cluster_name"
./${mode}-scaling.sh 1200 50 3 $ppn $nnodes weak $cluster_name