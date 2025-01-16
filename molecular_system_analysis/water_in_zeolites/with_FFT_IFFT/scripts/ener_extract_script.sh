#!/bin/bash

for i in {25..450..25};
do
	for j in RUN1 RUN2;
	do
		grep -A9 '<====  A V E R A G E S  ====>'  ../../../common/2PT/H2O/300K/$i/$j/prod1.log | tail -n1 | gawk '{print $2}' >> all_ener.txt
	done
	echo "exit $i"
done


	
