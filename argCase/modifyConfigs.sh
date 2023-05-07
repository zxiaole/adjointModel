#!/bin/bash
#
# Xiaole Zhang 2023.05.05
argNames=(strB blaTEMF ermB intl1 Staph sulII tetW tetG)
index=0
echo "Length:" ${#argNames[@]}
for ((i=0; i< ${#argNames[@]}-1; i++))
do
	echo ${argNames[$i]}
	echo ${argNames[$i+1]}
	grep ${argNames[$i]} -rl ./config
	sed -i 's/'${argNames[$i]}'/'${argNames[${i}+1]}'/gi'  `grep ${argNames[$i]} -rl ./config`
	rename 's/'${argNames[$i]}'/'${argNames[${i}+1]}'/' data/dep/*.bin
	polair3d config/caseConfig.cfg
	cp -a ./results/ ./results_${argNames[$i+1]}/
	rm ./results/*.*
done

#
#
