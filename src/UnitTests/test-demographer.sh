#!/bin/bash

epifile=../../inputs/inputs-test-demographer.txt
arvfile=../../inputs/inputs-art.txt

../ReferenceODE/simulate --year=2030 --epidemic=${epifile} --arvs=${arvfile} > output/ode-test-demographer.txt

i=0
while ((${i}<100))
do
  index=${i}
  while ((${#index}<2)); do index="0${index}"; done
  echo ${index}
  ./test-demographer --epidemic=${epifile} --years=52 > output/ibm-test-demographer-${index}.txt
  i=$((${i}+1))
done
