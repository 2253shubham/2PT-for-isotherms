#!/bin/bash

arr=()
while IFS= read -r line; do
    echo "Text read from file: $line"
    arr+=($(printf "%f" $line))
    #arr+=($line)
done <avg_all_ener.txt

cd ../output_files_2
declare -i flag=0 # default=0
for i in {25..450..25}; do
    echo "entering $i"
    cd $i
    python ../../code_files/IAG/therm2_prop_calc.py avg-output.txt ${arr[$flag]} $i ../100-Henry/avg_ref_dos.txt ../100-Henry/avg_ref_data.txt
    echo "executing code complete at this step"
    cd ../
    flag=$(($flag + 1))
    echo "exit $i"
done
