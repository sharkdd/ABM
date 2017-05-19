#!/bin/bash

epifile=../../inputs/inputs-test-progression.txt
arvfile=../../inputs/inputs-art.txt

./ode-test-progression --year=2030 --epidemic=${epifile} --arvs=${arvfile} > output/ode-test-progression.txt

i=0
while ((${i}<100))
do
  index=${i}
  while ((${#index}<2)); do index="0${index}"; done
  echo ${index}
  ./test-progression --epidemic=${epifile} --years=52 > output/ibm-test-progression-${index}.txt
  i=$((${i}+1))
done
