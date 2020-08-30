#!/bin/bash

if (( $# < 5 )); then
	echo Usage: $0 [ws_testbed] [dimensions] [ppn] [nnodes] [partition]
	exit 1
fi

ws_testbed=$1
dimensions=$2
ppn=$3
nnodes=$4
partition=$5
runs=1
heads=2
inst_directory="instances${dimensions}d"
outname="bm_inter_weak${dimensions}d_${partition}_$(date +%F_%H-%M)"

echo "Benchmark weak scaling internode"
echo "${dimensions}D version"
echo ""

echo "Instance,Nodes,NoShMem Efficiency,NoShMem GFlops,H=$heads s=X Efficiency,H=$heads s=X GFlops" > ${outname}.csv
for (( nodes=1; nodes<=$nnodes; nodes++ )) do
    p=$((nodes * ppn))
	instance=${ws_testbed}-$(printf "%02i" $nodes).txt
    printf "%s,%2i," $instance $nodes >> ${outname}.csv

    # NO SHARED MEMORY
    time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x		\
                --no-shared-memory                          \
                --input-instance=$inst_directory/$instance)
        timerstop=$(date +%s)
        echo "$output"
        echo "[*] mpirun took about $((timerstop - timerstart))s"
        echo ""
        perfresult=$(echo "$output" | grep Performance | grep -oe "[0-9]*\.[0-9]*")
        timeresult=$(echo "$output" | grep elapsed | grep -oe "[0-9]*\.[0-9]*")
        time=$(echo "scale=3; $time + $timeresult" | bc)
        perf=$(echo "scale=3; $perf + $perfresult" | bc)
    done

    # Computing efficiency
    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
    if [ $nodes -eq 1 ]; then
        baseline_noshm=$time
    fi
    efficiency=$(echo "scale=3; $baseline_noshm / $time" | bc)
	printf "%.3f," $efficiency >> ${outname}.csv
	printf "%.3f," $perf >> ${outname}.csv

    # SHARED MEMORY
    time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x				\
                        --split-direction=X                         \
                        --heads-per-shared-region=$heads			\
                        --input-instance=$inst_directory/$instance)
        timerstop=$(date +%s)
        echo "$output"
        echo "[*] mpirun took about $((timerstop - timerstart))s"
        echo ""
        perfresult=$(echo "$output" | grep Performance | grep -oe "[0-9]*\.[0-9]*")
        timeresult=$(echo "$output" | grep elapsed | grep -oe "[0-9]*\.[0-9]*")
        time=$(echo "scale=3; $time + $timeresult" | bc)
        perf=$(echo "scale=3; $perf + $perfresult" | bc)
    done

    # Computing efficiency
    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
    if [ $nodes -eq 1 ]; then
        baseline_shm=$time
    fi
    efficiency=$(echo "scale=3; $baseline_shm / $time" | bc)
	printf "%.3f," $efficiency >> ${outname}.csv
	printf "%.3f" $perf >> ${outname}.csv

    echo >> ${outname}.csv
done

sleep 5

printf "\n${dimensions}D Internode Weak scaling on $partition ($runs runs per method)\n" > ${outname}.table
column -s, -t < ${outname}.csv >> ${outname}.table
echo "END"
