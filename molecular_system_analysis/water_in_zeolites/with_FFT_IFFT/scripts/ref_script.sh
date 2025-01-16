#!/bin/bash

#!/bin/bash

arr=()
while IFS= read -r line; do
    echo "Text read from file: $line"
    arr+=($(printf "%f" $line))
    #arr+=($line)
done < ref_ener.txt

declare -i flag=0
echo "entering outside loop"
for j in RUN1 RUN2;
do
	echo $flag
	echo ${arr[$flag]}
	echo "linking start"
	cd ../2PT_H2O_files/100-Henry/$j
	ln -s /home/common/2PT/H2O/300K/100-Henry/$j/prod1* .
	ln -s /home/common/2PT/H2O/300K/100-Henry/analysis.tpr .
	echo "linking complete"
	cd ../../../without_FFT_IFFT/output_files/ref_data
	echo "executing code"
	python ../../dos_H2O_ref.py ../../../2PT_H2O_files/100-Henry/$j/analysis.tpr ../../../2PT_H2O_files/100-Henry/$j/prod1.trr ../../../2PT_H2O_files/100-Henry/$j/prod1.xtc 100-Henry-$j-tr_dos_ref.txt 100-Henry-$j-tr_ref_data.txt ${arr[$flag]} $j 
	echo "executing code complete at this step"
	cd ../../
    flag=$(($flag+1))
done
echo "exit outside the loop"

cd output_files/ref_data
python ../../avg.py 100-Henry-RUN1-tr_dos_ref.txt 100-Henry-RUN2-tr_dos_ref.txt avg_ref_tr_dos.txt
python ../../avg.py 100-Henry-RUN1-tr_ref_data.txt 100-Henry-RUN2-tr_ref_data.txt avg_ref_data.txt
python ../../avg.py tot_corr_per_10_RUN1.txt tot_corr_per_10_RUN2.txt avg_tot_corr_per_10.txt
python ../../avg.py dos_per_10_RUN1.txt dos_per_10_RUN2.txt avg_dos_per_10.txt

cd ../../