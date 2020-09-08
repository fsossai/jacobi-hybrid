#!/bin/bash

if (( $# < 6 )); then
    echo Usage: $0 <dimensions> <start_size> <iterations> <ppn> <nnodes> <partition>
    exit 1
fi

dimensions=$1
start_size=$2
iterations=$3
ppn=$4
nnodes=$5
partition=$6
runs=1
outname="bm_${partition}_topo${dimensions}d_$(date +%F_%H-%M)"
separator="|,"
other="Y"
if [ $dimensions -eq 3 ]; then
	other="Z"
fi

echo "Benchmarking the effects of different domain decompositions on $cluster_name"
echo "${dimensions}D version"
echo ""

echo "Instance,Processes,${separator}P:Topology,P:Performance,P:Time,${separator}H:Topology,H:Performance,H:Time" > ${outname}.csv
s_max=$(awk '{ print int( ($1 * $2) ^ (1 / $3) ) }' <<< "$nnodes $ppn $dimensions")
for (( s=1; s<=$s_max; s++ )) do
    p=$(awk '{ printf "%i", ($1 ^ $2) }' <<< "$s $dimensions")
    size=$(awk '{ printf "%.0f", $1 * ($2 ^ (1 / $3)) }' <<< "$start_size $p $dimensions" )

    # Logging problem size
    printf "%i" $size >> ${outname}.csv
    for (( i=1; i<$dimensions; i++ )) do
        printf "x%i" $size >> ${outname}.csv
    done
    printf ",%2i,${separator}" $p >> ${outname}.csv

    # Preparing instance
	instance="${size}\n"
    for (( i=1; i<$dimensions; i++ )) do
        instance="${instance}${size}\n"
    done
    instance="${instance}0.8\n1.0\n1e-9\n50\n"
	instance=$(echo -e "$instance")

    # Multidimensional decomposition
    time=0
	perf=0
	for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x		\
                --no-shared-memory							\
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

	# Writing results
    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
    topology=$(echo "$output" | grep topology | egrep -oe "[0-9]+(x[0-9]*)*")
	printf "%s," $topology >> ${outname}.csv
	printf "%.3f," $perf >> ${outname}.csv
	printf "%.3f,${separator}" $time >> ${outname}.csv

    # ------------------------------------------------------------------------

	# Unidimensional decomposition
	time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x         	\
                                --no-shared-memory					\
                                --cart-no-split-direction=Y			\
                                --cart-no-split-direction=$other	\
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

	# Writing results
    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
    topology=$(echo "$output" | grep topology | egrep -oe "[0-9]+(x[0-9]*)*")
	printf "%s," $topology >> ${outname}.csv
	printf "%.3f," $perf >> ${outname}.csv
	printf "%.3f\n" $time >> ${outname}.csv
done

sleep 5

printf "\n${dimensions}D Topology scaling on $partition ($runs runs per method)\n" > ${outname}.table
column -s, -t < ${outname}.csv >> ${outname}.table
echo "END"
