#!/bin/bash

if (( $# < 5 )); then
    echo "Usage: $0 [name] [dimension] [start_size] [nodes] [iterations]"
    exit 1
fi

name=$1
dimension=$2
start_size=$3
nodes=$4
iterations=$5

echo name=$name
echo dimension=$dimension
echo start_size=$start_size
echo nodes=$nodes
echo iterations=$iterations

for (( n=1; n<=$nodes; n++ )) do
    n_str=$(printf "%02i" $n)
    new_size=$(awk '{ printf("%.0f", $1 * ($2 ^ (1 / $3)) ) }' <<< "$start_size $n $dimension" )
    for (( i=1; i<=$dimension; i++ )) do
        echo $new_size >> ${name}-${n_str}.txt
    done
    echo 0.8 >> ${name}-${n_str}.txt
    echo 1.0 >> ${name}-${n_str}.txt
    echo 1e-9 >> ${name}-${n_str}.txt
    echo $iterations >> ${name}-${n_str}.txt
done

echo Done.

