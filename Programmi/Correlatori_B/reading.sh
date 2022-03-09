#!/bin/bash


awk '$9==15 || $9==1 {n=21}{if(n>0)n--; else if($1!="#")print $0}' /Users/manuel/Misure_Bz_Fixed/T20L24_beta3.787_b93/meson_corrs >> /Users/manuel/Misure_Bz_Fixed/T20L24_beta3.787_b93/meson_corrs_Cleaned_326
 
#awk '$9==15 || $9==1 {n=41}{if(n>0)n--; else if($1!="#")print $0}' meson_corrs_Bz93 >> /Users/manuel/Documents/GitHub/Conducibility/Programmi/Correlatori_B/meson_corrs_Bz93_input.out


