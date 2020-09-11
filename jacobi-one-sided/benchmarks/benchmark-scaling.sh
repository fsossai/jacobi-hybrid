#!/bin/bash

if (( $# < 5 )); then
    echo "Usage: $0 <cluster_name> <ppn> <nnodes> <mode> <scaling>"
    echo "mode: intranode | internode"
    echo "scaling: strong | weak"
    exit 1
fi

cluster_name=$1
ppn=$2
nnodes=$3
mode=$4
scaling=$5

if [[ ! $scaling =~ ^(strong|weak) ]]; then
    echo "ERROR: <scaling> must be 'strong' or 'weak'"
    exit 1
fi

# Set here the problem sizes

if [ "$mode" = "intranode" ]; then
    nnodes=""
    iterations=50
    if [ "$scaling" = "strong" ]; then
        start_size_2d=20000
        start_size_3d=700
    elif [ "$scaling" = "weak" ]; then
        start_size_2d=4000
        start_size_3d=250
    fi
elif [ "$mode" = "internode" ]; then
    nnodes=$SLURM_JOB_NUM_NODES
    iterations=50
    if [ "$scaling" = "strong" ]; then
        start_size_2d=44000
        start_size_3d=1200
    elif [ "$scaling" = "weak" ]; then
        start_size_2d=44000
        start_size_3d=1200
    fi
else
    echo "ERROR: Wrong argument for <mode> ($mode)"
    echo "mode: intranode | internode"
    exit 1
fi

echo "Running ${mode}-scaling.sh (2D / $scaling) on $cluster_name"
./${mode}-scaling.sh 2 $start_size_2d $iterations $ppn $nnodes $scaling $cluster_name

echo "Running ${mode}-scaling.sh (3D / $scaling) on $cluster_name"
./${mode}-scaling.sh 3 $start_size_3d $iterations $ppn $nnodes $scaling $cluster_name
