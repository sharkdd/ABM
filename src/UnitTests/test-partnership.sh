#!/bin/bash

inputs=../../inputs/inputs-test-partnership.txt

../pair-ode/pairs --config=${inputs} --year=2010 > output/ode-test-partnership.txt

i=0
while ((${i}<100))
do
  index=${i}
  while ((${#index}<2)); do index="0${index}"; done
  echo ${index}
  ./test-partnership --epidemic=${inputs} --years=32 > output/ibm-test-partnership-${index}.txt
  i=$((${i}+1))
done
