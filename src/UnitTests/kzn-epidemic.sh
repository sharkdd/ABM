#!/bin/bash

## Performs k repetitions using the IBM; defaults to 10 if no argument
## given.
k=${1:-10}

i=0
while ((${i}<${k}))
do

  index="${i}"
  while ((${#index}<${#k})); do index="0${index}"; done

  echo "Running IBM replication ${index}"
  ../simulate --art --epidemic=../../inputs/test-kzn-epidemic.txt --analysis=../../inputs/inputs-art.txt --years=52 > output/ibm-kzn-epi-${index}.txt

  i=$((${i}+1))

done
echo "Running reference ODE simulation"
../ReferenceODE/simulate --art --epidemic=../../inputs/test-kzn-epidemic.txt --arvs=../../inputs/inputs-art.txt --year=2030 > output/ode-kzn-epi.txt
