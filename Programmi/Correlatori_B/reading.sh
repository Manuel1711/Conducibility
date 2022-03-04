#!/bin/bash


awk '$9==15 || $9==1 {n=41}{if(n>0)n--; else if($1!="#")print $0}' /Users/manuel/Misure_Bz_Fixed/T40L48_beta4.14_b0/meson_corrs >> /Users/manuel/Misure_Bz_Fixed/T40L48_beta4.14_b0/meson_corrs_Cleaned_47
 
#awk '$9==15 || $9==1 {n=41}{if(n>0)n--; else if($1!="#")print $0}' meson_corrs_Bz93 >> /Users/manuel/Documents/GitHub/Conducibility/Programmi/Correlatori_B/meson_corrs_Bz93_input.out


