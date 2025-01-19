#!/bin/bash

arr=()
while IFS= read -r line; do
    echo "Text read from file: $line"
    arr+=($(printf "%f" $line))
    #arr+=($line)
done <all_ener.txt

declare -i flag=28
for i in {375..450..25}; do # modify it based on systems to investigate
    echo "entering $i"
    for j in RUN1 RUN2; do # modify it based on number of runs
        echo $flag
        echo ${arr[$flag]}
        echo "linking start"
        cd ../2PT_H2O_files/$i/$j                         # check locations
        ln -s /home/common/2PT/H2O/300K/$i/$j/prod1* .    # check locations
        ln -s /home/common/2PT/H2O/300K/$i/analysis.tpr . # check locations
        echo "linking complete"
        MEMORY=$(ps -eo pid,ppid,cmd,%mem,%cpu --sort=-%mem | head | awk '{print $0}')
        echo "$MEMORY"
        cd ../../../with_FFT_IFFT/output_files/$i
        echo "executing code"
        python ../../dos_H2O.py ../../../2PT_H2O_files/$i/$j/analysis.tpr ../../../2PT_H2O_files/$i/$j/prod1.trr ../../../2PT_H2O_files/$i/$j/prod1.xtc $j-output.txt ${arr[$flag]} $j ../ref_data/avg_ref_tr_dos.txt ../ref_data/avg_ref_data.txt
        #python ../../dos_H2O_2.py ../../2PT_H2O_files/$i/$j/analysis.tpr ../../2PT_H2O_files/$i/$j/prod1.trr ../../2PT_H2O_files/$i/$j/prod1.xtc $j-output.txt ${arr[$flag]} $j ../ref_data/avg_ref_tr_dos.txt ../ref_data/avg_ref_data.txt
        echo "executing code complete at this step"
        cd ../../
        flag=$(($flag + 1))
    done
    echo "exit $i"
    MEMORY = $(ps -eo pid,ppid,cmd,%mem,%cpu --sort=-%mem | head | awk '{print $0}')
    echo "$MEMORY"

    cd output_files/$i
    python ../../avg.py tot_corr_per_10_for_RUN1.txt tot_corr_per_10_for_RUN2.txt avg_tot_corr_per_10.txt
    python ../../avg.py dos_per_10_for_RUN1.txt dos_per_10_for_RUN2.txt avg_dos_per_10.txt
    cd ../../

done
