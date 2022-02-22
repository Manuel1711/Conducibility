#!/bin/bash

set -e

PATH=$PATH:/Users/manuel/Documents/GitHub/Conducibility/Programmi/

optimization="0 -g"



g++ -O$optimization -std=c++11 -o Diff_Means Diff_Means.C `root-config --libs --cflags` -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include -L/usr/local/include  -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas

./Diff_Means
