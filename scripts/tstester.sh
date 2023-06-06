#!/bin/bash

# "1e-9" "1e-8"  TODO
dts=("1e-7" "1e-6" "1e-5" "1e-4" "1e-3" "1e-2" "3e-7" "3e-6" "3e-5" "3e-4" \
     "3e-3" "3e-2" "1.0" "0.3" "1e-1" "3e-8" "1e-8")

for dt in ${dts[*]}
do 
    echo "Starting run $dt"
    echo $dt > errorpropagation/results_$dt.csv
    time ./bin/tslosserrs $dt >> errorpropagation/results_$dt.csv
done
