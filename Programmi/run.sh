#!/bin/bash

set -e

PATH=$PATH:/Users/manuel/Documents/GitHub/Conducibility/Programmi/

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@             SMEARING EXTRACTION            @@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


#@@@ COMPLILAZIONE Spectral_gsl.C @@@#

optimization="0 -g"

#boolEsBin=0 non binnato (il valore di Es si sceglie in pars.h), boolEsBin=1 binnato 
boolEsBin=0

if [ $boolEsBin -eq 0 ]
then
    N="-D SINGLE"
fi

if [ $boolEsBin -eq 1 ]
then
    N="-D CICLO"
fi

if [ $boolEsBin -eq 2 ]
then
    N="-D CICLO_LAMBDA"
fi



#method
boolM=0
#basis
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



g++ -O$optimization -std=c++11 -o Smearing_Func Smearing_Func.C -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include -L/usr/local/include  -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas $Method $Basis $N




if [ $boolEsBin -eq 0 ]
then
    
    ./Smearing_Func
    
fi




if [ $boolEsBin -eq 1 ]
then

    
    
    #check e cancellazione file
    
    if test -f "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/bin.out";then
	rm /Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/bin.out
    fi
    if test -f "/tmp/temp.out";then
	rm /tmp/temp.out
    fi
    if test -f "/tmp/rho_o.out";then
	rm /tmp/rho_o.out
    fi
    if test -f "/tmp/rho_e.out";then
	rm /tmp/rho_e.out
    fi

    #apro i file per le funzioni spettrali
    
    if [ $boolEO -eq 0 ]
    then
	echo "@type xydy" >> /tmp/rho_e.out
    fi
    if [ $boolEO -eq 1 ]
    then
	echo "@type xydy" >> /tmp/rho_o.out
    fi

    #inizio ciclo sugli Estar
    
    for k in `seq 0 1`
    do
	
	for i in `seq 0 9`
	do
	    if [ $k -eq 0 ]
	    then
		b="0.$i"
	    fi
	    if [ $k -eq 1 ]
	    then
		b="1.$i"
	    fi
	    
	    echo $b >> /Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/bin.out
	    
	    ./Smearing_Func
	    
	    rm /Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/bin.out

	    #leggo rho ed errore associato
	    read r err < /tmp/temp.out

	    #scrivo sul file aperto insieme ai tempi euclidei
	    if [ $boolEO -eq 0 ]
	    then
		echo $b $r $err >> /tmp/rho_e.out
	    fi
	    if [ $boolEO -eq 1 ]
	    then
		echo $b $r $err >> /tmp/rho_o.out
	    fi
	    
	    rm /tmp/temp.out
	    
	    
	done
    done
    
fi


