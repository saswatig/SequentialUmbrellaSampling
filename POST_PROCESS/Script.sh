#!/bin/bash
echo ../RESULTS/Histo_UM > INPUT
echo LOG_PCHI > OUTPUT
ifort -O3 UM_LN_HISTO.f
./a.out

