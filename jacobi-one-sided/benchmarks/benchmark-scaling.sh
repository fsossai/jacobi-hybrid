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
    start_size_2d=9000
    start_size_3d=400
elif [ "$mode" = "internode" ]; then
	nnodes=$SLURM_JOB_NUM_NODES
    start_size_2d=44000
    start_size_3d=1200
else
    echo "ERROR: Wrong argument for <mode> ($mode)"
    echo "mode: intranode | internode"
    exit 1
fi

echo "Running ${mode}-scaling.sh (2D / strong) on $cluster_name"
./${mode}-scaling.sh 2 $start_size_2d 50 $ppn $nnodes strong $cluster_name

echo "Running ${mode}-scaling.sh (3D / strong) on $cluster_name"
./${mode}-scaling.sh 3 $start_size_3d 50 $ppn $nnodes strong $cluster_name

echo "Running ${mode}-scaling.sh (2D / weak) on $cluster_name"
./${mode}-scaling.sh 2 $start_size_2d 50 $ppn $nnodes weak $cluster_name

echo "Running ${mode}-scaling.sh (3D / weak) on $cluster_name"
./${mode}-scaling.sh 3 $start_size_3d 50 $ppn $nnodes weak $cluster_name