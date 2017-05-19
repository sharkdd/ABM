#!/bin/bash

epifile=../../inputs/inputs-test-transmission.txt

./test-transmission --epidemic=${epifile} > output/ibm-test-transmission.txt
