#!/bin/bash

if (( $# < 3 )); then
	echo "Usage: $0 <cluster_name> <ppn> <nnodes>"
	exit 1
fi

cluster_name=$1
ppn=$2
nnodes=$3

echo "Running topology-scaling.sh 2D on $cluster_name"
./topology-scaling.sh 2 44000 50 $ppn $nnodes $cluster_name

echo "Running topology-scaling.sh 3D on $cluster_name"
./topology-scaling.sh 3 1200 50 $ppn $nnodes $cluster_name