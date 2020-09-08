#!/bin/bash

if (( $# < 2 )); then
	echo "Usage: $0 <cluster_name> <max_procs>"
	exit 1
fi

cluster_name=$1
max_procs=$2

echo "Running topology-scaling.sh 2D on $cluster_name"
./topology-scaling.sh 2 44000 50 $max_procs $cluster_name

echo "Running topology-scaling.sh 3D on $cluster_name"
./topology-scaling.sh 3 1200 50 $max_procs $cluster_name