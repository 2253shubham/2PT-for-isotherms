#!/bin/bash

# modify the code for appropriate file locations

#for i in 100-Henry;
for i in {25..450..25}; do
    echo "entering $i"
    for j in RUN1 RUN2; do
        echo "executing code"
        cd output_files/$i
        python ../../dos_computation.py ../../../2PT_H2O_files/$i/$j/analysis.tpr ../../../2PT_H2O_files/$i/$j/prod1.trr ../../../2PT_H2O_files/$i/$j/prod1.xtc $j
        cd ../../
        #flag=$(($flag+1))
    done
    cd output_files/ref_data
    python ../../avg.py tot_dos_per_10_for_RUN1.txt tot_dos_per_10_for_RUN2.txt avg_tot_dos_per_10_for-100-Henry-system.txt
    python ../../avg.py tr_dos_per_10_for_RUN1.txt tr_dos_per_10_for_RUN2.txt avg_tr_dos_per_10_for-100-Henry-system.txt
    python ../../avg.py rot_dos_per_10_for_RUN1.txt rot_dos_per_10_for_RUN2.txt avg_rot_dos_per_10_for-100-Henry-system.txt
    cd ../../
done
