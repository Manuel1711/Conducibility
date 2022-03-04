#include <iostream>
#include <math.h>
#include "mp.h"
#include "statistical.h"
#include <random>

using namespace std;

int NMass=3, Dim=3, Nt=40,NBz=3;
int NConf[3]={47,102,56};
int Nboot=100;

string M[3] = {"mu", "md", "ms"};
string D_s[3] = {"xx", "yy", "zz"};
string Dir_z[2] = {"xy", "z"};
int Bz[3] = {0,41,93};


Real C[3] = {0.6666666, -0.33333333, -0.33333333};

Real e2 = 4*Pi/(137.04);


int main(){


  Real Corr_bz[3][NConf[1]][NMass][Dim][Nt];

  Real *tE = new Real[Nt];
  
  FILE *Correlators_bz_Inputs;
  char open_Correlators_bz_Inputs[1024];
  
  
  
  for(int bz=0; bz<NBz; bz++){
    
    sprintf(open_Correlators_bz_Inputs, "/Users/manuel/Misure_Bz/T40L48_beta4.14_b%d/meson_corrs_Cleaned_%d", Bz[bz], NConf[bz]);
    if ((Correlators_bz_Inputs = fopen(open_Correlators_bz_Inputs, "r")) == NULL ){
      printf("Error opening the input file: %s\n",open_Correlators_bz_Inputs);
      exit(EXIT_FAILURE);
    }
    for(int iconf=0; iconf<NConf[bz]; iconf++){
      for(int d=0; d<Dim; d++){
	for(int m=0; m<NMass; m++){
	  for(int t=0; t<Nt; t++){
	    
	    char trash1[1024], trash2[1024];
	    double trash3;
	    
	    fscanf(Correlators_bz_Inputs, "%s " "%s " "%lf\n", trash1, trash2, &trash3);
	    Corr_bz[bz][iconf][m][d][t] = conv(trash2);
	    if(bz==0){
	      tE[t] = conv(trash1);
	      //cout << "t: " << tE[t] << " C: " << Corr_bz[bz][iconf][m][d][t] << endl;
	    }
	    
	    
	  }//t
        }//d
      }//m
    }//iconf
    
  }//bz
  
  fclose(Correlators_bz_Inputs);

  
  Real Corr_bz_m[3][NConf[1]][Dim][Nt];
  for(int bz=0; bz<NBz; bz++){
    for(int iconf=0; iconf<NConf[bz]; iconf++){
      for(int d=0; d<Dim; d++){
	for(int t=0; t<Nt; t++){

	  for(int m=0; m<NMass; m++){
	    
	    Corr_bz_m[bz][iconf][d][t] += e2*C[m]*C[m]/16*Corr_bz[bz][iconf][m][d][t];
	    
	  }
	  
	}
      }
    }
  }
  
  Real Corr_bz_m_d[3][NConf[1]][2][Nt];
  for(int bz=0; bz<NBz; bz++){
    for(int iconf=0; iconf<NConf[bz]; iconf++){
      for(int t=0; t<Nt; t++){
	
	for(int d=0; d<Dim; d++){
	  
	  if(bz==0) Corr_bz_m_d[bz][iconf][0][t] += Corr_bz_m[bz][iconf][d][t];
	  else if(bz>0){
	    if(d==0 or d==1){ 
	      Corr_bz_m_d[bz][iconf][0][t] += Corr_bz_m[bz][iconf][d][t];
	    }
	    else if(d==2) Corr_bz_m_d[bz][iconf][1][t] = Corr_bz_m[bz][iconf][d][t];
	  }
	}
	if(bz==0) Corr_bz_m_d[bz][iconf][0][t] = Corr_bz_m_d[bz][iconf][0][t]/3;  
	else if(bz>0) Corr_bz_m_d[bz][iconf][0][t] = Corr_bz_m_d[bz][iconf][0][t]/2;
	
      }
    }
  }

  
  Real Corr_bz_m_d_mu[3][2][Nt];   

  for(int bz=0; bz<NBz; bz++){
    for(int d=0; d<2; d++){
      for(int t=0; t<Nt; t++){
	
	for(int iconf=0; iconf<NConf[bz]; iconf++){
	  Corr_bz_m_d_mu[bz][d][t] += Corr_bz_m_d[bz][iconf][d][t];
	}

	Corr_bz_m_d_mu[bz][d][t] = Corr_bz_m_d_mu[bz][d][t]/NConf[bz];
	cout << "bz: " << bz << " d " << d << " t: " << t << " C: " << Corr_bz_m_d_mu[bz][d][t] << endl;
      }
    }
  }
  
  
  Real Corr_bz_m_d_err[3][2][Nt];
  
  for(int bz=0; bz<NBz; bz++){
    for(int d=0; d<2; d++){
      for(int t=0; t<Nt; t++){
	
	for(int iconf=0; iconf<NConf[bz]; iconf++){
	  Corr_bz_m_d_err[bz][d][t] += pow(Corr_bz_m_d[bz][iconf][d][t] - Corr_bz_m_d_mu[bz][d][t],2);
	}
	Corr_bz_m_d_err[bz][d][t] = sqrt(Corr_bz_m_d_err[bz][d][t])/NConf[bz];
      }
    }
  }
  
  
 Real Diff[3][NConf[1]][2][Nt];

 for(int bz=0; bz<NBz; bz++){
   for(int d=0; d<2; d++){
     for(int iconf=0; iconf<NConf[bz]; iconf++){
       for(int t=0; t<Nt/2-1; t++){
	 Diff[bz][iconf][d][t]= Corr_bz_m_d[bz][iconf][d][t+1] - Corr_bz_m_d[bz][iconf][d][39-t];
	 cout << "KKK: " << Diff[bz][iconf][d][t] << "  " << t+1 << "  " << 39-t << endl;
       }
     }
   }
 }


 Real Diff_M[3][2][Nt/2-1];
 
 for(int bz=0; bz<NBz; bz++){
   for(int d=0; d<2; d++){
     for(int t=0; t<Nt/2-1; t++){

       for(int iconf=0; iconf<NConf[bz]; iconf++){
	 Diff_M[bz][d][t] += Diff[bz][iconf][d][t];
       }
       Diff_M[bz][d][t] = Diff_M[bz][d][t]/NConf[bz];
     }
   }
 }

 Real Diff_Err[3][2][Nt/2-1];

 for(int bz=0; bz<NBz; bz++){
   for(int d=0; d<2; d++){
     for(int t=0; t<Nt/2-1; t++){

       for(int iconf=0; iconf<NConf[bz]; iconf++){
	 Diff_Err[bz][d][t] += pow(Diff[bz][iconf][d][t] - Diff_M[bz][d][t],2);
       }
       
       Diff_Err[bz][d][t] = sqrt(Diff_Err[bz][d][t])/NConf[bz];
       
     }
   }
 }
 
 FILE *Correlators_bz_outputs_M, *Correlators_bz_outputs_Means;
 char open_Correlators_bz_outputs_M[1024], open_Correlators_bz_outputs_Means[1024];
 
 for(int bz=0; bz<NBz; bz++){
   for(int d=0; d<2; d++){
     sprintf(open_Correlators_bz_outputs_M, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Differenze/T40L48_beta4.14_b%d/Diff_b%d_M_%s.out", Bz[bz], Bz[bz], Dir_z[d].c_str());
     if ((Correlators_bz_outputs_M = fopen(open_Correlators_bz_outputs_M, "w")) == NULL ){
       printf("Error opening the input file: %s\n",open_Correlators_bz_outputs_M);
       exit(EXIT_FAILURE);
     }
     sprintf(open_Correlators_bz_outputs_Means, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Differenze/T40L48_beta4.14_b%d/Means_b%d_M_%s.out", Bz[bz], Bz[bz], Dir_z[d].c_str());
     if ((Correlators_bz_outputs_Means = fopen(open_Correlators_bz_outputs_Means, "w")) == NULL ){
       printf("Error opening the input file: %s\n",open_Correlators_bz_outputs_Means);
       exit(EXIT_FAILURE);
     }
     
     for(int t=0; t<Nt; t++){
       
       if(t<Nt/2-1) fprintf(Correlators_bz_outputs_M, "%d " "%s " " %s\n", t, conv(abs(Diff_M[bz][d][t])).c_str(), conv(Diff_Err[bz][d][t]).c_str());
       fprintf(Correlators_bz_outputs_Means, "%d " "%s " " %s\n", t, conv(Corr_bz_m_d_mu[bz][d][t]).c_str(), conv(Corr_bz_m_d_err[bz][d][t]).c_str());
       
       
     }
   }
   
 }

 fclose(Correlators_bz_outputs_M);
 
 return 0;
 
}
