#!/bin/bash

arr=()
while IFS= read -r line; do
    echo "Text read from file: $line"
    arr+=($(printf "%f" $line))
    #arr+=($line)
done <../avg_all_ener.txt

cd ../../output_files_hs
declare -i flag=0 # default=0
for i in {25..450..25}; do
    echo "entering $i"
    cd $i
    python ../../code_files/hs/therm2_prop_calc.py avg-output-hs.txt ${arr[$flag]} $i
    echo "executing code complete at this step"
    cd ../
    flag=$(($flag + 1))
    echo "exit $i"
done
