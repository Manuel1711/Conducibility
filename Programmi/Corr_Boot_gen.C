#include <iostream>
#include <math.h>
#include "mp.h"
#include <random>

using namespace std;

int NMass=3, Dim=3, Nt=40,NBz=3;

int Nboot = 100;
int NConf[3]={47,101,56};

string M[3] = {"mu", "md", "ms"};
string D_s[3] = {"xx", "yy", "zz"};
string Dir_z[2] = {"xy", "z"};
int Bz[3] = {0,41,93};


Real C[3] = {0.6666666, -0.33333333, -0.33333333};

Real e2 = 4*Pi/(137.04);


int main(){
  

  //Input correlators for each value of the magnetic field
  
  Real *****Corr_bz = new Real****[3];
  for(int bz=0; bz<3; bz++){
    Corr_bz[bz] = new Real***[NConf[1]];
    for(int iconf = 0; iconf < 101; iconf++){
      Corr_bz[bz][iconf] = new Real**[NMass];
      for(int m = 0; m < NMass; m++){
	Corr_bz[bz][iconf][m] = new Real*[Dim];
	for(int d=0; d<Dim; d++){
	  Corr_bz[bz][iconf][m][d] = new Real[Nt];
	}
      }
    }
  }
  

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
      for(int m=0; m<NMass; m++){
	for(int d=0; d<Dim; d++){
	  for(int t=0; t<Nt; t++){
	    
	    char trash1[1024], trash2[1024];
	    double trash3;
	     
	    fscanf(Correlators_bz_Inputs, "%s " "%s " "%lf\n", trash1, trash2, &trash3);
	    Corr_bz[bz][iconf][m][d][t] = conv(trash2);
	    if(bz==0){
	      tE[t] = conv(trash1);
	      cout << "t: " << tE[t] << " C: " << Corr_bz[bz][iconf][m][d][t] << endl;
	    }
	    
	    
	  }//t
        }//d
      }//m
    }//iconf
    
  }//bz
  
  fclose(Correlators_bz_Inputs);

  
  //Generation of bootstrap events
  Real ******Corr_bz_boot = new Real*****[3];
  for(int bz=0; bz<3; bz++){
  Corr_bz_boot[bz] = new Real****[NConf[1]];
    for(int iconf = 0; iconf < NConf[1]; iconf++){
      Corr_bz_boot[bz][iconf] = new Real***[NMass];
      for(int m = 0; m < NMass; m++){
	Corr_bz_boot[bz][iconf][m] = new Real**[Dim];
	for(int d=0; d<Dim; d++){
	  Corr_bz_boot[bz][iconf][m][d] = new Real*[Nt];
	  for(int t=0; t<Nt; t++){
	    Corr_bz_boot[bz][iconf][m][d][t] = new Real[Nboot];
	  }
	}
      }
    }
  }
  random_device rd;
  mt19937 gen(rd());
  for(int bz=0; bz<NBz; bz++){
    uniform_int_distribution<> distrib(0,NConf[bz]-1);
    for(int m=0; m<NMass; m++){
      for(int d=0; d<Dim; d++){
	for(int t=0; t<Nt; t++){
	  for(int iboot=0; iboot<Nboot; iboot++){
	    for(int iconf=0; iconf<NConf[bz]; iconf++){
	      //cout << "INT: " << distrib(gen) << endl; 
	      Corr_bz_boot[bz][iconf][m][d][t][iboot] = Corr_bz[bz][distrib(gen)][m][d][t];
	    } 
	  }
	  
	}
      }
    }
  }
  for(int bz=0; bz<3; bz++){
    for(int iconf=0; iconf<4; iconf++){
      for(int m = 0; m < NMass; m++){
	for(int d=0; d<Dim; d++){
	  delete[] Corr_bz[bz][iconf][m][d];
	}
	delete[] Corr_bz[bz][iconf][m];
      }
      delete[] Corr_bz[bz][iconf];
    }
    delete[] Corr_bz[bz];
  }
  delete[] Corr_bz;


  
  //Media sulle configurazioni
  Real *****Corr_bz_boot_m = new Real****[3];
  for(int bz=0; bz<3; bz++){
    Corr_bz_boot_m[bz] = new Real***[NMass];
    for(int m = 0; m < NMass; m++){
      Corr_bz_boot_m[bz][m] = new Real**[Dim];
      for(int d=0; d<Dim; d++){
	Corr_bz_boot_m[bz][m][d] = new Real*[Nt];
	for(int t=0; t<Nt; t++){
	  Corr_bz_boot_m[bz][m][d][t] = new Real[Nboot];
	}
      }
    }
  }
  for(int bz=0; bz<3; bz++){
    for(int m = 0; m < NMass; m++){
      for(int d=0; d<Dim; d++){
	for(int t=0; t<Nt; t++){
	  for(int iboot=0; iboot<Nboot; iboot++){
	    for(int iconf=0; iconf<NConf[bz]; iconf++){
	      Corr_bz_boot_m[bz][m][d][t][iboot] += Corr_bz_boot[bz][iconf][m][d][t][iboot];
	    }
	    Corr_bz_boot_m[bz][m][d][t][iboot]=Corr_bz_boot_m[bz][m][d][t][iboot]/NConf[bz];
	  }
	}
      }
    }
  }
  for(int bz=0; bz<3; bz++){
    for(int iconf=0; iconf<4; iconf++){
      for(int m = 0; m < NMass; m++){
	for(int d=0; d<Dim; d++){
	  for(int t=0; t<Nt; t++){
	    delete[] Corr_bz_boot[bz][iconf][m][d][t];
	  }
	  delete[] Corr_bz_boot[bz][iconf][m][d];
	}
	delete[] Corr_bz_boot[bz][iconf][m];
      }
    delete[] Corr_bz_boot[bz][iconf];
    }
    delete[] Corr_bz_boot[bz];
  }
  delete[] Corr_bz_boot;


  
  //Somma sui 3 flavor pesato dalle cariche
  Real ****Corr_bz_boot_mf = new Real***[3];
  for(int bz=0; bz<3; bz++){
    Corr_bz_boot_mf[bz] = new Real**[Dim];
    for(int d=0; d<Dim; d++){
      Corr_bz_boot_mf[bz][d] = new Real*[Nt];
      for(int t=0; t<Nt; t++){
	Corr_bz_boot_mf[bz][d][t] = new Real[Nboot];
      }
    }
  }
  for(int bz=0; bz<3; bz++){
    for(int d=0; d<Dim; d++){
      for(int t=0; t<Nt; t++){
	for(int iboot=0; iboot<Nboot; iboot++){
	  
	  for(int m = 0; m < NMass; m++){
	    Corr_bz_boot_mf[bz][d][t][iboot] += C[m]*C[m]*Corr_bz_boot_m[bz][m][d][t][iboot]/16;
	  }
	  //cout << "bz: " << bz << " d: " << d << " t: " << t << " iboot: " << iboot << " C: " << Corr_bz_boot_mf[bz][d][t][iboot]<< endl;
	}
      }
    }
  }
  for(int bz=0; bz<3; bz++){
    for(int m = 0; m < NMass; m++){
      for(int d=0; d<Dim; d++){
	for(int t=0; t<Nt; t++){
	  delete[] Corr_bz_boot_m[bz][m][d][t];
	}
	delete[] Corr_bz_boot_m[bz][m][d];
      }
      delete[] Corr_bz_boot_m[bz][m];
    }
    delete[] Corr_bz_boot_m[bz];
  }
  delete[] Corr_bz_boot_m;
  

  
  //Media su x,y
  Real ****Corr_bz_boot_mf_xz = new Real***[3];
  for(int bz=0; bz<3; bz++){
    Corr_bz_boot_mf_xz[bz] = new Real**[Dim-1];
    for(int d=0; d<Dim-1; d++){
      Corr_bz_boot_mf_xz[bz][d] = new Real*[Nt];
      for(int t=0; t<Nt; t++){
	Corr_bz_boot_mf_xz[bz][d][t] = new Real[Nboot];
      }
    }
  }
  for(int bz=0; bz<3; bz++){
    for(int t=0; t<Nt; t++){
      for(int iboot=0; iboot<Nboot; iboot++){
	
	for(int d=0; d<Dim; d++){
	  if(d==0 or d==1){
	    Corr_bz_boot_mf_xz[bz][0][t][iboot] += Corr_bz_boot_mf[bz][d][t][iboot];
	  }
	  else if(d==2) Corr_bz_boot_mf_xz[bz][1][t][iboot] = Corr_bz_boot_mf[bz][d][t][iboot];
        }
	cout << "BBBB: " << Corr_bz_boot_mf_xz[bz][0][t][iboot] << endl;
	Corr_bz_boot_mf_xz[bz][0][t][iboot] = Corr_bz_boot_mf_xz[bz][0][t][iboot]/2;
      }
    }
  }
  for(int bz=0; bz<3; bz++){
    for(int d=0; d<Dim; d++){
      for(int t=0; t<Nt; t++){
	delete[] Corr_bz_boot_mf[bz][d][t];
      }
      delete[] Corr_bz_boot_mf[bz][d];
    }
    delete[] Corr_bz_boot_mf[bz];
  }
  delete[] Corr_bz_boot_mf;

  
  //Output bootstrap events
  FILE *Boot_output;
  char open_Boot_output[1024];
  
  for(int bz=0; bz<3; bz++){
    for(int d=0; d<Dim-1; d++){
      
      sprintf(open_Boot_output,"/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Not_Subtracted/Bootstraps/T20L24_beta3.787_b%d/boot_correlators_%s.out", Bz[bz], Dir_z[d].c_str());
      if ((Boot_output = fopen(open_Boot_output, "w")) == NULL ){
	printf("Error opening the input file: %s\n",open_Boot_output);
	exit(EXIT_FAILURE);
      }
      for(int t=0; t<Nt; t++){
	for(int iboot=0; iboot<Nboot; iboot++){

	  fprintf(Boot_output, "%s " "%s\n", conv(tE[t]).c_str(), conv(Corr_bz_boot_mf_xz[bz][d][t][iboot]).c_str());
	  
	}
      }
    }
  }

  fclose(Boot_output);
  
  
  //Computing the related \mu and \sigma bootstrap
  FILE *MuSigma_output;
  char open_MuSigma_output[1024];
  
  for(int bz=0; bz<3; bz++){
    for(int d=0; d<Dim-1; d++){
      
      sprintf(open_MuSigma_output,"/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Not_Subtracted/Means_Sigmas/T20L24_beta3.787_b%d/mu_sigma_correlators_%s.out", Bz[bz], Dir_z[d].c_str());
      if ((MuSigma_output = fopen(open_MuSigma_output, "w")) == NULL ){
	printf("Error opening the input file: %s\n",open_MuSigma_output);
	exit(EXIT_FAILURE);
      }
      
      for(int t=0; t<Nt; t++){
	Real Mu=0;
	Real Bar_Sigma=0;
	for(int iboot=0; iboot<Nboot; iboot++){
	  Mu += Corr_bz_boot_mf_xz[bz][d][t][iboot];
	  cout << "bz: " << bz << " d: " << d << " t: " << t << " iboot: " << iboot << " Mu: " << Corr_bz_boot_mf_xz[bz][d][t][iboot] << endl;
	  Bar_Sigma += Corr_bz_boot_mf_xz[bz][d][t][iboot]*Corr_bz_boot_mf_xz[bz][d][t][iboot];
	}

	fprintf(MuSigma_output, "%s " "%s " "%s\n", conv(tE[t]).c_str(), conv(Mu/Nboot).c_str(), conv((Bar_Sigma-Mu*Mu/Nboot)/Nboot).c_str());
	
      }
    }
  }

  fclose(MuSigma_output);
  


  
  //Subtraction of bz=0
  Real ****Corr_bz_boot_mf_xz_sub = new Real***[2];
  for(int bz=0; bz<2; bz++){
  Corr_bz_boot_mf_xz_sub[bz] = new Real**[Dim-1];
    for(int d=0; d<Dim-1; d++){
      Corr_bz_boot_mf_xz_sub[bz][d] = new Real*[Nt];
      for(int t=0; t<Nt; t++){
	Corr_bz_boot_mf_xz_sub[bz][d][t] = new Real[Nboot];
      }
    }
  }
  for(int bz=0; bz<2; bz++){
    for(int d=0; d<Dim-1; d++){
      for(int t=0; t<Nt; t++){
	for(int iboot=0; iboot<Nboot; iboot++){

	  Corr_bz_boot_mf_xz_sub[bz][d][t][iboot] = Corr_bz_boot_mf_xz[bz+1][d][t][iboot]-Corr_bz_boot_mf_xz[0][d][t][iboot];
	  
	}
      }
    }
  }
  for(int bz=0; bz<2; bz++){
    for(int d=0; d<Dim-1; d++){
      for(int t=0; t<Nt; t++){
	delete[] Corr_bz_boot_mf_xz[bz][d][t];
      }
      delete[] Corr_bz_boot_mf_xz[bz][d];
    }
    delete[] Corr_bz_boot_mf_xz[bz];
  }
  delete[] Corr_bz_boot_mf_xz;
  
  
  //Output bootstrap events subtracted correlators 
  FILE *Boot_sub_output;
  char open_Boot_sub_output[1024];
  
  for(int bz=0; bz<2; bz++){
    for(int d=0; d<Dim-1; d++){
      
      sprintf(open_Boot_sub_output,"/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Bootstraps/T20L24_beta3.787_b%d_Subtracted/boot_correlators_%s.out", Bz[bz+1], Dir_z[d].c_str());
      if ((Boot_sub_output = fopen(open_Boot_sub_output, "w")) == NULL ){
	printf("Error opening the input file: %s\n",open_Boot_sub_output);
	exit(EXIT_FAILURE);
      }
      for(int t=0; t<Nt; t++){
	for(int iboot=0; iboot<Nboot; iboot++){
	  
	  fprintf(Boot_sub_output, "%s " "%s\n", conv(tE[t]).c_str(), conv(Corr_bz_boot_mf_xz_sub[bz][d][t][iboot]).c_str());
	  
	}
      }
    }
  }

  fclose(Boot_sub_output);
  
  
  //Computing the related subtracted \mu and \sigma bootstrap
  FILE *MuSigma_sub_output;
  char open_MuSigma_sub_output[1024];
  
  for(int bz=0; bz<2; bz++){
    for(int d=0; d<Dim-1; d++){
      
      sprintf(open_MuSigma_sub_output,"/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Means_Sigmas/T20L24_beta3.787_b%d_Subtracted/mu_sigma_correlators_%s.out", Bz[bz+1], Dir_z[d].c_str());
      if ((MuSigma_sub_output = fopen(open_MuSigma_sub_output, "w")) == NULL ){
	printf("Error opening the input file: %s\n",open_MuSigma_sub_output);
	exit(EXIT_FAILURE);
      }
      
      for(int t=0; t<Nt; t++){
	Real Mu=0;
	Real Bar_Sigma=0;
	for(int iboot=0; iboot<Nboot; iboot++){
	  Mu += Corr_bz_boot_mf_xz_sub[bz][d][t][iboot];
	  Bar_Sigma += Corr_bz_boot_mf_xz_sub[bz][d][t][iboot]*Corr_bz_boot_mf_xz_sub[bz][d][t][iboot];
	}
	
	fprintf(MuSigma_sub_output, "%s " "%s " "%s\n", conv(tE[t]).c_str(), conv(Mu/Nboot).c_str(), conv((Bar_Sigma-Mu*Mu/Nboot)/Nboot).c_str());
	
      }
    }
  }

  fclose(MuSigma_sub_output);
  
  
  return 0;
  
}
