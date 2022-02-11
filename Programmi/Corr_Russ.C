#include <iostream>
#include <math.h>
#include "mp.h"
 
using namespace std;

#define NMass 3
#define Dim 3
#define Nt 40

int NConf=108;

string M[3] = {"mu", "md", "ms"};
string D_s[3] = {"xx", "yy", "zz"};

Real C[3] = {0.6666666, -0.33333333, -0.33333333};

Real e2 = 4*Pi/(137.04);

int main(){

  Real Corr[NConf][NMass][Dim][Nt];
  Real tE[Nt];

  FILE *Correlators_Inputs;
  char open_Correlators_Inputs[1024];
  
  
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Misure/L48_T12_b4.06100000_ml0.0009168830_ms0.0258102551_zero_strange/nissa_mu0.200/meson_corrs_Cleaned");
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Misure/L48_T10_b3.85_ml0.0014_ms0.0394_zero_strange/nissa_mu0.200/meson_corrs_Cleaned");
  sprintf(open_Correlators_Inputs, "/Users/manuel/Misure/mu_zero/L48_T10_b3.96127000_ml0.0011192722_ms0.0315075218/meson_corrs_Cleaned");
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Misure/L48_T12_b3.93998000_ml0.0011677267_ms0.0328715394_zero_strange/nissa_mu0.140/meson_corrs_Cleaned");
  
  if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs);
    exit(EXIT_FAILURE);
  }
  
  char trash1[1024], trash2[1024];
  double trash3;
  
  for(int conf=0; conf<NConf; conf++){ 
    for(int m=0; m<NMass; m++){
      for(int d=0; d<Dim; d++){
	
	cout << "conf: " << conf << " d: " << d << " m: " << m << endl;
	for(int t=0; t<Nt; t++){
	  
	    fscanf(Correlators_Inputs, "%s " "%s " "%lf\n", trash1, trash2, &trash3);
	    tE[t] = conv(trash1);
	    Corr[conf][m][d][t] = conv(trash2);
	    cout << "tE: " << tE[t] << " Corr: " << Corr[conf][m][d][t] << endl;
	    
	    
	}//t
  
      }//d
    }//m
  }//conf
  
  
  //Media
  Real CorrMu[NMass][Dim][Nt];
  for(int m=0; m<NMass; m++){
    for(int d=0; d<Dim; d++){
      for(int t=0; t<Nt; t++) {
	for(int conf=0; conf<NConf; conf++){
	  
	  CorrMu[m][d][t] += Corr[conf][m][d][t];
	  
	}
	
	CorrMu[m][d][t] = CorrMu[m][d][t]/NConf;
	
      }
    }
  }
  
  //Sigma
  
  Real CorrErr[NMass][Dim][Nt];
  for(int m=0; m<NMass; m++){
    for(int d=0; d<Dim; d++){
      for(int t=0; t<Nt; t++){
	for(int conf=0; conf<NConf; conf++){
	  
	  CorrErr[m][d][t] += pow(Corr[conf][m][d][t] - CorrMu[m][d][t], 2);
	  
	}
	CorrErr[m][d][t] = sqrt(CorrErr[m][d][t]/pow(NConf,2));
      }
    }
  }
  

  //Output Correlators Jackknife
  FILE *Corr_Jack_Out;
  char open_Corr_Jack_Out[1024];
  
  //Sum of the correlators weighted by their charges
  Real CorrMu_Tot[Dim][Nt];
  Real CorrErr_Tot[Dim][Nt];

  for(int d=0; d<Dim; d++){
    for(int t=0; t<Nt; t++) {
      
      for(int m=0; m<NMass; m++){
	
	CorrMu_Tot[d][t] += C[m]*C[m]*CorrMu[m][d][t]/16;
	CorrErr_Tot[d][t] = pow(C[m]*C[m],2)/16*CorrErr[m][d][t]*CorrMu[m][d][t];
	
      }
      CorrErr_Tot[d][t] = sqrt(CorrErr_Tot[d][t]);
      
    }
  }

  
  
  //Sum of the correlators on the three directions
  Real CorrMu_TotMean[Nt], CorrErr_TotMean[Nt];
  for(int t=0; t<Nt; t++){
    
    for(int d=0; d<Dim; d++){  
      
      CorrMu_TotMean[t] += CorrMu_Tot[d][t];
      CorrErr_TotMean[t] += CorrErr_Tot[d][t]*CorrErr_Tot[d][t];
    }
    
    
    CorrMu_TotMean[t] =  4*e2*CorrMu_TotMean[t];
    CorrErr_TotMean[t] = 4*e2*sqrt(CorrErr_TotMean[t]);
    cout << t << "  " << CorrMu_TotMean[t] <<  "  " << CorrErr_TotMean[t] << endl;

    //C'è un fattore 12 di differenza
    
  }

  
  
  
  return 0;
  
}
