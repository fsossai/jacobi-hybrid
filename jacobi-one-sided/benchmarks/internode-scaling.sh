#!/bin/bash

if (( $# < 7 )); then
	echo "Usage: $0 <dimensions> <start_size> <iterations> <ppn> <nnodes> <mode> <partition>"
    echo "mode: strong | weak"
	exit 1
fi

dimensions=$1
start_size=$2
iterations=$3
ppn=$4
nnodes=$5
mode=$6
partition=$7
runs=1
heads=2
separator="|," 

if [[ ! $mode =~ ^(strong|weak) ]]; then
    echo "ERROR: <mode> must be 'strong' or 'weak'"
    exit 1
fi
outname="bm_${partition}_inter_${mode}${dimensions}d_$(date +%F_%H-%M)"

echo "Benchmark: ${mode} scaling internode"
echo "${dimensions}D version, using $heads heads per node"
echo ""

echo "Instance,Nodes,${separator}P:Topology,P:Efficiency,P:Performance,P:Time,${separator}H:Topology,H:Efficiency,H:Performance,H:Time" > ${outname}.csv
for (( nodes=1; nodes<=$nnodes; nodes++ )) do
    p=$((nodes * ppn))

    if [ "$mode" = "strong" ]; then
        size=$start_size
    elif [ "$mode" = "weak" ]; then
	    size=$(awk '{ printf "%.0f", $1 * ($2 ^ (1 / $3)) }' <<< "$start_size $nodes $dimensions" )
    fi

	# Logging problem size and nodes
    printf "(%i)^%i," $size $dimensions >> ${outname}.csv
    printf "%i,${separator}" $nodes >> ${outname}.csv

    # Preparing instance
    instance="${size}\n"
    for (( i=1; i<$dimensions; i++ )) do
        instance="${instance}${size}\n"
    done
    instance="${instance}0.8\n1.0\n1e-9\n${iterations}\n"
    instance=$(echo -e "$instance")

    # NO SHARED MEMORY
    time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x		\
                --no-shared-memory                          \
                <<< "$instance")
        timerstop=$(date +%s)
        echo "$output"
        echo "[*] mpirun took about $((timerstop - timerstart))s"
        echo ""
        perfresult=$(echo "$output" | grep Performance | grep -oe "[0-9]*\.[0-9]*")
        timeresult=$(echo "$output" | grep elapsed | grep -oe "[0-9]*\.[0-9]*")
        time=$(echo "scale=3; $time + $timeresult" | bc)
        perf=$(echo "scale=3; $perf + $perfresult" | bc)
    done

    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
	topology=$(echo "$output" | grep topology | egrep -oe "[0-9]+(x[0-9]*)*")
    if [ $nodes -eq 1 ]; then
        baseline_pure=$time
    fi

    # Computing efficiency
    if [ "$mode" = "strong" ]; then
        efficiency=$(echo "scale=3; $baseline_pure / $time / $nodes" | bc)
    elif [ "$mode" = "weak" ]; then
        efficiency=$(echo "scale=3; $baseline_pure / $time" | bc)
    fi

	# Writing results
	printf "%s," $topology >> ${outname}.csv
	printf "%.3f," $efficiency >> ${outname}.csv
	printf "%.3f," $perf >> ${outname}.csv
	printf "%.3f,${separator}" $time >> ${outname}.csv

	# ---------------------------------------------------------------------------------

    # SHARED MEMORY
    time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x				\
                        --split-direction=X                         \
                        --heads-per-shared-region=$heads			\
                        <<< "$instance")
        timerstop=$(date +%s)
        echo "$output"
        echo "[*] mpirun took about $((timerstop - timerstart))s"
        echo ""
        perfresult=$(echo "$output" | grep Performance | grep -oe "[0-9]*\.[0-9]*")
        timeresult=$(echo "$output" | grep elapsed | grep -oe "[0-9]*\.[0-9]*")
        time=$(echo "scale=3; $time + $timeresult" | bc)
        perf=$(echo "scale=3; $perf + $perfresult" | bc)
    done

    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
	topology=$(echo "$output" | grep topology | egrep -oe "[0-9]+(x[0-9]*)*")
    if [ $nodes -eq 1 ]; then
        baseline_hybrid=$time
    fi

    # Computing efficiency
    if [ "$mode" = "strong" ]; then
        efficiency=$(echo "scale=3; $baseline_hybrid / $time / $nodes" | bc)
    elif [ "$mode" = "weak" ]; then
        efficiency=$(echo "scale=3; $baseline_hybrid / $time" | bc)
    fi

	# Writing results
	printf "%s," $topology >> ${outname}.csv
	printf "%.3f," $efficiency >> ${outname}.csv
	printf "%.3f," $perf >> ${outname}.csv
	printf "%.3f" $time >> ${outname}.csv

    echo >> ${outname}.csv
done

sleep 5

printf "\n${dimensions}D Internode $mode scaling on $partition ($runs runs per method)\n" > ${outname}.table
column -s, -t < ${outname}.csv >> ${outname}.table
echo "END"
