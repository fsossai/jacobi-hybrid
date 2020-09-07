#!/bin/bash

if (( $# < 6 )); then
	echo Usage: $0 [ws_testbed] [dimensions] [ppn] [hpn] [nnodes] [partition]
	exit 1
fi

ws_testbed=$1
dimensions=$2
ppn=$3
hpn=$4
nnodes=$5
partition=$6
runs=1
inst_directory="instances${dimensions}d"
outname="bm_topo_weak${dimensions}d_${partition}_$(date +%F_%H-%M)"
other="Y"
if [ $dimensions -eq 3 ]; then
	other="Z"
fi

echo "Benchmark Weak scaling internode with squared cartesian topologies"
echo "${dimensions}D version"
echo ""

echo "Instance,Nodes,Multidim Topology,Multidim Performance,Unidim Topology,Unidim Performance" > ${outname}.csv
s_max=$(awk '{ printf("%i", int( ($1*$2) ^ (1 / $3) / 2 )) }' <<< "$nnodes $hpn $dimensions")
for (( s=1; s<=$s_max; s++ )) do
	nodes=$(echo "scale=0; (($hpn * $s) ^ $dimensions) / 2" | bc)
    p=$(echo "scale=0; $ppn * $nodes" | bc)
	instance=${ws_testbed}-$(printf "%02i" $nodes).txt
    printf "%s,%3i," $instance $nodes >> ${outname}.csv

	# Multidimensional decomposition
    time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x		\
				--heads-per-shared-region=$hpn				\
                --input-instance=$inst_directory/$instance)
        timerstop=$(date +%s)
        echo "$output"
        echo "[*] mpirun took about $((timerstop - timerstart))s"
        echo ""
        perfresult=$(echo "$output" | grep Performance | grep -oe "[0-9]*\.[0-9]*")
        timeresult=$(echo "$output" | grep elapsed | grep -oe "[0-9]*\.[0-9]*")
		topology=$(echo "$output" | grep topology | egrep -oe "[0-9]+(x[0-9]*)*")
        time=$(echo "scale=3; $time + $timeresult" | bc)
        perf=$(echo "scale=3; $perf + $perfresult" | bc)
    done

	# Writing results
	printf "%s," $topology >> ${outname}.csv
	printf "%.3f," $perf >> ${outname}.csv

	# Unidimensional decomposition
    time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x		\
				--cart-no-split-direction=Y					\
				--cart-no-split-direction=$other			\
				--heads-per-shared-region=$hpn				\
                --input-instance=$inst_directory/$instance)
        timerstop=$(date +%s)
        echo "$output"
        echo "[*] mpirun took about $((timerstop - timerstart))s"
        echo ""
        perfresult=$(echo "$output" | grep Performance | grep -oe "[0-9]*\.[0-9]*")
        timeresult=$(echo "$output" | grep elapsed | grep -oe "[0-9]*\.[0-9]*")
		topology=$(echo "$output" | grep topology | egrep -oe "[0-9]+(x[0-9]*)*")
        time=$(echo "scale=3; $time + $timeresult" | bc)
        perf=$(echo "scale=3; $perf + $perfresult" | bc)
    done

	# Writing results
	printf "%s," $topology >> ${outname}.csv
	printf "%.3f\n" $perf >> ${outname}.csv
done

sleep 5

printf "\n${dimensions}D Internode Weak scaling with squared cartesian topologies on $partition ($runs runs per method)\n" > ${outname}.table
column -s, -t < ${outname}.csv >> ${outname}.table
echo "END"
