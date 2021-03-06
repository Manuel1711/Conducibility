#!/bin/bash

set -e

PATH=$PATH:/Users/manuel/Documents/GitHub/Conducibility/Programmi/


optimization="0 -g"

boolM=1
boolB=1
#declared="-D ABAB"

if [ $boolM -eq 0 ]
then
    Method="-D BG"
fi

if [ $boolM -eq 1 ]
then
    Method="-D HLN"
fi


if [ $boolB -eq 0 ]
then
    Basis="-D EXP"
fi

if [ $boolB -eq 1 ]
then
    Basis="-D COS"
fi

g++ -O$optimization -std=c++11 -o Corr_Recon Corr_Recon.C -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include -L/usr/local/include  -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas $Method $Basis


./Corr_Recon
