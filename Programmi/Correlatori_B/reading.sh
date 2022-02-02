#!/bin/bash


awk '$9==15 || $9==1 {n=11}{if(n>0)n--; else if($1!="#")print $0}' /Users/manuel/Misure/mu_zero/L48_T10_b3.96127000_ml0.0011192722_ms0.0315075218/meson_corrs >> /Users/manuel/Misure/mu_zero/L48_T10_b3.96127000_ml0.0011192722_ms0.0315075218/meson_corrs_Cleaned

#awk '$9==15 || $9==1 {n=41}{if(n>0)n--; else if($1!="#")print $0}' meson_corrs_Bz93 >> /Users/manuel/Documents/GitHub/Conducibility/Programmi/Correlatori_B/meson_corrs_Bz93_input.out


