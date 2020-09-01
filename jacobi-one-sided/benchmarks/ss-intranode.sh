#!/bin/bash

if (( $# < 4 )); then
	echo Usage: $0 [instance] [dimensions] [ppn] [partition]
	exit 1
fi

instance=$1
dimensions=$2
ppn=$3
partition=$4
runs=1
inst_directory="instances${dimensions}d"
outname="bm_intra_strong${dimensions}d_${partition}_$(date +%F_%H-%M)"

echo "Instance,Processes,NoShMem Speedup,NoShMem GFlops,ShMem s=X Speedup,ShMem s=X GFlops" > ${outname}.csv
for (( p=1; p<=$ppn; p++ )) do
    printf "%s,%2i," $instance $p >> ${outname}.csv
    
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

    # Computing speedup
    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
    if [ $p -eq 1 ]; then
        baseline_noshm=$time
    fi
    speedup=$(echo "scale=3; $baseline_noshm / $time" | bc)
    echo -n $speedup "," >> ${outname}.csv
	printf "%.3f," $perf >> ${outname}.csv

    # SHARED MEMORY
    time=0
	perf=0
    for (( i=1; i<=$runs; i++ )) do
        timerstart=$(date +%s)
        output=$(mpirun -n $p ../jacobi${dimensions}d.x				\
                        --split-direction=X                         \
                        --heads-per-shared-region=1					\
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

    # Computing speedup
    time=$(echo "scale=3; $time / $runs" | bc)
    perf=$(echo "scale=3; ($perf / $runs) / 1000" | bc)
    if [ $p -eq 1 ]; then
        baseline_shm=$time
    fi
    speedup=$(echo "scale=3; $baseline_shm / $time" | bc)
    echo -n $speedup "," >> ${outname}.csv
	printf "%.3f" $perf >> ${outname}.csv

    echo >> ${outname}.csv
done

sleep 5

printf "\n${dimensions}D Intranode Strong scaling on $partition ($runs runs per method)\n" > ${outname}.table
column -s, -t < ${outname}.csv >> ${outname}.table
echo "END"
