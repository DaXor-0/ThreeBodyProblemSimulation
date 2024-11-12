#!/bin/bash


# Check if enough arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 N WARMUP NOFBODIES ITERATIONS NOFPROC"
    echo "N: Total number of runs"
    echo "WARMUP: Number of warmup runs to ignore in statistics"
    exit 1
fi

# Read input arguments
N=$1
WARMUP=$2
NOFBODIES=$3
ITERATIONS=$4
NOFPROC=$5

mpi='yes'
mpi='no'


# Check that N is greater than WARMUP
if [ "$N" -le "$WARMUP" ]; then
    echo "Error: Total runs (N) must be greater than warm-up runs (WARMUP)"
    exit 1
fi

# Array to store execution times
times=()

# Function to calculate mean
calculate_mean() {
    sum=0
    for time in "${times[@]}"; do
        sum=$(echo "$sum + $time" | bc -l)
    done
    echo "scale=5; $sum / ${#times[@]}" | bc -l
}

# Function to calculate variance and standard deviation
calculate_variance_stddev() {
    mean=$1
    sum_sq_diff=0
    for time in "${times[@]}"; do
        diff=$(echo "$time - $mean" | bc -l)
        sq_diff=$(echo "$diff * $diff" | bc -l)
        sum_sq_diff=$(echo "$sum_sq_diff + $sq_diff" | bc -l)
    done
    variance=$(echo "scale=5; $sum_sq_diff / ${#times[@]}" | bc -l)
    stddev=$(echo "scale=5; sqrt($variance)" | bc -l)
    echo "$variance $stddev"
}

# Run the program N times
for (( i=1; i<=N; i++ )); do
    start_time=$(date +%s.%N)
    
    if [ $mpi == 'yes' ]; then
        srun -n $NOFPROC ./simulation.out $NOFBODIES $ITERATIONS benchmark.csv
    else
        ./simulation.out $NOFBODIES $ITERATIONS benchmark.csv
    fi

    end_time=$(date +%s.%N)

    # Calculate elapsed time
    elapsed=$(echo "$end_time - $start_time" | bc -l)

    # Print elapsed time for each run
    echo "Run $i: $elapsed seconds"

    # Only store times after warmup period
    if [ "$i" -gt "$WARMUP" ]; then
        times+=("$elapsed")
    fi
done

# Calculate average time
mean=$(calculate_mean)

# Calculate variance and standard deviation
read variance stddev <<< $(calculate_variance_stddev "$mean")

# Display the results
echo "Average time (after warmup): $mean seconds"
echo "Variance: $variance"
echo "Standard Deviation: $stddev"
