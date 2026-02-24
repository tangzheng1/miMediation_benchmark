#!/bin/bash

# Get 8 parameters
template=$1
n=$2
p=$3
num1_A=$4
num1_B=$5
num2=$6
d=$7
seed=$8

# Print basic info
echo "=========================================="
echo "Template: $template | n=$n p=$p | A=$num1_A B=$num1_B num2=$num2 | d=$d seed=$seed"
echo "=========================================="


# Run R script
Rscript global_test.R "$template" "$n" "$p" "$num1_A" "$num1_B" "$num2" "$d" "$seed"

exit_code=$?

# Report result
if [ $exit_code -eq 0 ]; then
    echo "✓ Success (seed=$seed)"
else
    echo "✗ Failed with exit code $exit_code"
fi

exit $exit_code
