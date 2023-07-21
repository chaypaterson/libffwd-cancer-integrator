#!/bin/bash

# compute the root mean squared error from a time series of local root squared
# error (std. deviations)

for file in results_*; 
do 
    step=$(sed '1q' $file); 
    err=$(awk 'BEGIN {FS = ","}; 
               (NR > 2){totsqe += $2 * $2; n += 1} 
               END {print sqrt(totsqe / n)}' $file);
    echo "$step, $err,"
done
