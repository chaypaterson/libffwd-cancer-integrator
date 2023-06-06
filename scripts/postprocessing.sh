#!/bin/bash

for file in results_*; 
do 
    step=$(sed '1q' $file); 
    err=$(awk '(max < $2 && NR > 2){max = $2} END {print max}' $file);
    echo "$step, $err,"
done
