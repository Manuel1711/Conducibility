#include <iostream>
#include <math.h>
#include "mp.h"

using namespace std;

#define NMass 3
#define Dim 3
#define Nt 40

int NConf[2] = {102,56};

//#define NConf 102
int Bz[2] = {41,93};

string M[3] = {"mu", "md", "ms"};
string D_s[3] = {"xx", "yy", "zz"};
Real C[3] = {0.6666666, -0.33333333, -0.33333333};
//#conf per ogni cluster
int nc = 12;
//#cluster
int NClust[2]={10,5};
//Nc = 10 per Bz=102 (9 gruppi da 10(nc=10) e 1 da 12(nc=12)), Nc=5 per Bz=56 (4 gruppi da 11(nc=11), 1 da 12(nc=12)).


int main(){
  
  
  Real Corr[2][NConf[1]][NMass][Dim][Nt];
  Real tE[Nt];

  
  
  FILE *Correlators_Inputs;
  char open_Correlators_Inputs[1024];

  
  for(int B = 0; B<2; B++){
    
    sprintf(open_Correlators_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Correlatori_B/meson_corrs_Bz%d_input.out", Bz[B]);
    
    
    if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
      printf("Error opening the input file: %s\n",open_Correlators_Inputs);
      exit(EXIT_FAILURE);
    }
    
    char trash1[1024], trash2[1024];
    double trash3;
    for(int conf=0; conf<NConf[B]; conf++){ 
      for(int d=0; d<Dim; d++){
	for(int m=0; m<NMass; m++){
	  cout << "conf: " << conf << " d: " << d << " m: " << m << endl;
	  for(int t=0; t<Nt; t++){
	    
	    fscanf(Correlators_Inputs, "%s " "%s " "%lf\n", trash1, trash2, &trash3);
	    tE[t] = conv(trash1);
	    Corr[B][conf][m][d][t] = conv(trash2);
	    cout << "tE: " << tE[t] << " Corr: " << Corr[B][conf][m][d][t] << endl;
	    
	    
	  }//t
	}//d
      }//m
    }//conf
  }//B
  fclose(Correlators_Inputs);
  
    
  //Media
  Real CorrMu[2][NMass][Dim][Nt];
  for(int B = 0; B<2; B++){
    for(int m=0; m<NMass; m++){
      for(int d=0; d<Dim; d++){
	for(int t=0; t<Nt; t++) {
	  for(int conf=0; conf<NConf[B]; conf++){
	    
	    CorrMu[B][m][d][t] += Corr[B][conf][m][d][t];
	    
	  }
	
	  CorrMu[B][m][d][t] = CorrMu[B][m][d][t]/NConf[B];
	  
	}
      }
    }
  }//B
  
  //Sigma
  
  Real CorrErr[2][NMass][Dim][Nt];
  for(int B = 0; B<2; B++){
    for(int m=0; m<NMass; m++){
      for(int d=0; d<Dim; d++){
	for(int t=0; t<Nt; t++){
	  for(int conf=0; conf<NConf[B]; conf++){
	  
	    CorrErr[B][m][d][t] += pow(Corr[B][conf][m][d][t] - CorrMu[B][m][d][t], 2);
	    
	  }
	  CorrErr[B][m][d][t] = sqrt(CorrErr[B][m][d][t]/pow(NConf[B],2));
	}
      }
    }
  }
  //Output Correlatori

  FILE *Correlators_Outputs;
  char open_Correlators_Outputs[1024];

  for(int B = 0; B<2; B++){

    for(int m=0; m<NMass; m++){
      
      sprintf(open_Correlators_Outputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/Corr_%s_Bz%d.out", M[m].c_str(), Bz[B]);
      
      if ((Correlators_Outputs = fopen(open_Correlators_Outputs, "w")) == NULL ){
	printf("Error opening the input file: %s\n",open_Correlators_Outputs);
	exit(EXIT_FAILURE);
      }
      
      fprintf(Correlators_Outputs,"@type xydy");
    
      for(int d=0; d<Dim; d++){
	
	fprintf(Correlators_Outputs, "\n \n ");
	
	for(int t=0; t<Nt; t++){
	  
	  fprintf(Correlators_Outputs, "%s " "%s " "%s\n", conv(tE[t]).c_str(), conv(CorrMu[B][m][d][t]).c_str(), conv(CorrErr[B][m][d][t]).c_str());
	  
	  
	}
	
      }
      
    }
  }

  fclose(Correlators_Outputs);

  
  //Sum of the correlators weighted by their charges
  Real Corr_Tot[2][NConf[1]][Dim][Nt];
  
  for(int B = 0; B<2; B++){
    for(int conf=0; conf<NConf[B]; conf++){
      for(int d=0; d<Dim; d++){
	for(int t=0; t<Nt; t++) {

	  for(int m=0; m<NMass; m++){
	    
	    Corr_Tot[B][conf][d][t] += C[m]*Corr[B][conf][m][d][t];
	  } 
	  
	}
      }
    }
  }
  
  
  //Jackknife Sampling

  //Mean for each sample 
  Real Clust_Mu_41[NClust[0]][Dim][Nt], Clust_Mu_93[NClust[1]][Dim][Nt];
  for(int B = 0; B<2; B++){
    for(int d=0; d<Dim; d++){
      for(int t=0; t<Nt; t++){
	
	int inc=0, np=0, iconf=0;
	for(int i=0; i<NClust[B]; i++){
	  
	  if(i<NClust[B]-1){
	    if(B==0) inc=10;
	    if(B==1) inc=11;
	  }
	  else if(i==NClust[B]-1) inc = nc;
	  
	  for(int j=0; j<inc; j++){
	    iconf=j+np;
	    if(B==0)Clust_Mu_41[i][d][t] += Corr_Tot[B][iconf][d][t]; 
	    if(B==1)Clust_Mu_93[i][d][t] += Corr_Tot[B][iconf][d][t];
	  }
	  if(B==0)Clust_Mu_41[i][d][t]=Clust_Mu_41[i][d][t]/inc;
	  if(B==1)Clust_Mu_93[i][d][t]=Clust_Mu_93[i][d][t]/inc;
	  np=inc;
	}
	
      }
    }
  }
  
  
  //Jackknife averages
  Real Corr_Mu_Jack_41[NClust[0]][Dim][Nt], Corr_Mu_Jack_93[NClust[1]][Dim][Nt];
  
  for(int B = 0; B<2; B++){
    for(int d=0; d<Dim; d++){
      
      
      for(int t=0; t<Nt; t++){

	
	for(int i=0; i<NClust[B];i++){
	  
	  for(int j=0; j<NClust[B]; j++){
	    
	    
	    if(B==0)
	      if(j != i) Corr_Mu_Jack_41[i][d][t] += Clust_Mu_41[j][d][t];
	    if(B==1)
	      if(j != i) Corr_Mu_Jack_93[i][d][t] += Clust_Mu_93[j][d][t];
	    
	  }
	  if(B==0) Corr_Mu_Jack_41[i][d][t]=Corr_Mu_Jack_41[i][d][t]/(NClust[B]-1);
	  if(B==1) Corr_Mu_Jack_93[i][d][t]= Corr_Mu_Jack_93[i][d][t]/(NClust[B]-1);
	}
	
      }
    }
  }
  
  
  //Output Correlators Jackknife
  FILE *Corr_Jack_Out;
  char open_Corr_Jack_Out[1024];
  
  for(int B = 0; B<2; B++){
    for(int d=0; d<Dim; d++){
      
      sprintf(open_Corr_Jack_Out, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/Jack_Corr/Corr_Jack_%d_%s", Bz[B], D_s[d].c_str());
      
      if ((Corr_Jack_Out = fopen(open_Corr_Jack_Out, "w")) == NULL ){
	printf("Error opening the output file: %s\n",open_Corr_Jack_Out);
	exit(EXIT_FAILURE);
      }

      for(int i=0; i<NClust[B];i++){
	
	for(int t=0; t<Nt; t++){
	  
	  if(B==0) fprintf(Corr_Jack_Out, "%s "  "%s\n", conv(tE[t]).c_str(), conv(Corr_Mu_Jack_41[i][d][t]).c_str());
	  if(B==1) fprintf(Corr_Jack_Out, "%s "  "%s\n", conv(tE[t]).c_str(), conv(Corr_Mu_Jack_93[i][d][t]).c_str());
	  
	}

	fprintf(Corr_Jack_Out, " \n" );

      }
	
    }
  }
  
  
  //Jackknife sigmas
  Real Sigma_41[Dim][Nt], Sigma_93[Dim][Nt];
  for(int B = 0; B<2; B++){
    for(int d=0; d<Dim; d++){
      for(int t=0; t<Nt; t++){
	
	Real App=0, App_Sq=0;
	for(int i=0; i<NClust[B]; i++){

	  if(B==0){
	    App += Corr_Mu_Jack_41[i][d][t];
	    App_Sq += pow(Corr_Mu_Jack_41[i][d][t],2);
	  }
	  if(B==1){
	    App += Corr_Mu_Jack_93[i][d][t];
	    App_Sq += pow(Corr_Mu_Jack_93[i][d][t],2);
	  }

	}
	
	App = App/NClust[B];
	App_Sq = App_Sq/NClust[B];   

	if(B==0)Sigma_41[d][t] = sqrt(NClust[B]-1)*sqrt(App_Sq - pow(App,2));
	if(B==1)Sigma_93[d][t] = sqrt(NClust[B]-1)*sqrt(App_Sq - pow(App,2));
	
      }
    }
  }
  

  //Sigmas Output
  
  FILE *Sigma_Jack_Out;
  char open_Sigma_Jack_Out[1024];

  for(int B = 0; B<2; B++){
    for(int d=0; d<Dim; d++){
  
      sprintf(open_Sigma_Jack_Out, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/Jack_Corr/Sigmas_%d_%s", Bz[B], D_s[d].c_str());

      if ((Sigma_Jack_Out = fopen(open_Sigma_Jack_Out, "w")) == NULL ){
	printf("Error opening the output file: %s\n",open_Sigma_Jack_Out);
	exit(EXIT_FAILURE);
      }
      
      for(int t=0; t<Nt; t++){
	if(B==0)fprintf(Sigma_Jack_Out, "%s " "%s\n",conv(tE[t]).c_str(), conv(Sigma_41[d][t]).c_str());
	if(B==1)fprintf(Sigma_Jack_Out, "%s " "%s\n",conv(tE[t]).c_str(), conv(Sigma_93[d][t]).c_str());
      }

      fprintf(Sigma_Jack_Out, "\n");
      
    }

  }
      
  fclose(Sigma_Jack_Out);
  
  return 0;
  
}

