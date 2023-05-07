#!/bin/bash


echo ../LOG_PCHI > INPUT
echo ../../RESULTS/Normalized_PChi > OUTPUT
ifort -O3 NORM.f
./a.out

