#include <iostream>
#include <math.h>
#include "mp.h"
#include "statistical.h"
#include <random>


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
  
  Real *****Corr_bz = new Real****[3];
  for(int bz=0; bz<3; bz++){
    Corr_bz[bz] = new Real***[NConf[1]];
    for(int iconf = 0; iconf < NConf[1]; iconf++){
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
  
  Real *****Diff = new Real****[3];
  for(int bz=0; bz<3; bz++){
    Diff[bz] = new Real***[NConf[1]];
    for(int iconf = 0; iconf < NConf[1]; iconf++){
      Diff[bz][iconf] = new Real**[NMass];
      for(int m = 0; m < NMass; m++){
        Diff[bz][iconf][m] = new Real*[Dim];
        for(int d=0; d<Dim; d++){
          Diff[bz][iconf][m][d] = new Real[Nt/2-1];
        }
      }
    }
  }
  
  for(int bz=0; bz<NBz; bz++){
    for(int iconf=0; iconf<NConf[bz]; iconf++){
      for(int d=0; d<Dim; d++){
	for(int m=0; m<NMass; m++){
	  for(int t=0; t<Nt/2-1; t++){
	    
	    Diff[bz][iconf][m][d][t] = Corr_bz[bz][iconf][m][d][t+1]-Corr_bz[bz][iconf][m][d][39-t];
	    //if(bz==0 and d==0 and m==0 and t==0){
	      cout << t+1 << "   " << 39-t << endl;
	      cout << "AAAAA: " << Diff[bz][iconf][m][d][t] <<  "   "  << Corr_bz[bz][iconf][m][d][t+1] << "  " << Corr_bz[bz][iconf][m][d][39-t] <<  endl;
	      //}
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
  
  Real Media[3][Dim][NMass][Nt/2-1], Err[3][Dim][NMass][Nt/2-1];
  for(int bz=0; bz<NBz; bz++){
    for(int d=0; d<Dim; d++){
      for(int m=0; m<NMass; m++){
	for(int t=0; t<Nt/2-1; t++){ 
	  for(int iconf=0; iconf<NConf[bz]; iconf++){
	    
	    Media[bz][d][m][t]+= Diff[bz][iconf][m][d][t];
	    
	  }
	  cout << " PPP: " <<  Media[bz][d][m][t] << endl;
	  Media[bz][d][m][t] = Media[bz][d][m][t]/NConf[bz];
	  cout << " PPP2 : " <<  Media[bz][d][m][t] << endl;
	}
      }
    }
  }

  for(int bz=0; bz<NBz; bz++){
    for(int d=0; d<Dim; d++){
      for(int m=0; m<NMass; m++){
	for(int t=0; t<Nt/2-1; t++){
	  for(int iconf=0; iconf<NConf[bz]; iconf++){
	    
	    Err[bz][d][m][t] += pow(Diff[bz][iconf][m][d][t]-Media[bz][d][m][t],2);
	    //cout << "Diff: " << Diff[bz][iconf][m][d][t] << " M: " << Media[bz][d][m][t] << "  " << Err[bz][d][m][t] << endl;
	    cout  << Diff[bz][iconf][m][d][t]  << endl;
	  }

	  cout << endl;
	  Err[bz][d][m][t] = sqrt(Err[bz][d][m][t])/NConf[bz];
	  
	}
      }
    }
  }

  
  FILE *Correlators_bz_outputs;
  char open_Correlators_bz_outputs[1024];
  
  
  for(int bz=0; bz<NBz; bz++){

    sprintf(open_Correlators_bz_outputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Differenze/T40L48_beta4.14_b%d/Diff_b%d.out", Bz[bz], Bz[bz]);
    if ((Correlators_bz_outputs = fopen(open_Correlators_bz_outputs, "w")) == NULL ){
      printf("Error opening the input file: %s\n",open_Correlators_bz_outputs);
      exit(EXIT_FAILURE);
    }
    
    for(int d=0; d<Dim; d++){
      for(int m=0; m<NMass; m++){

	fprintf(Correlators_bz_outputs, "bz=%d " "d=%s " "m=%s\n", Bz[bz], D_s[d].c_str(), M[m].c_str());
	
	for(int t=0; t<Nt/2-1; t++){

	  fprintf(Correlators_bz_outputs,"%d " "%s " "%s\n", t, conv(Media[bz][d][m][t]).c_str(), conv(Err[bz][d][m][t]).c_str());
	  //cout << "bz: " << bz << " d: " << d << " m: " << m << " t: " << t << endl;
	  //cout << "D: " << Media[bz][d][m][t] << "  " << Err[bz][d][m][t] << endl;
	  
	}
	fprintf(Correlators_bz_outputs, "\n");
      }
    }
  }
  
  for(int bz=0; bz<3; bz++){
     for(int iconf=0; iconf<NConf[bz]; iconf++){
       for(int m = 0; m < NMass; m++){
	 for(int d=0; d<Dim; d++){
	   delete[] Diff[bz][iconf][m][d];
	 }
	 delete[]Diff[bz][iconf][m];
       }
       delete[] Diff[bz][iconf];
     }
     delete[] Diff[bz];
   }
  delete[] Diff;
  
  
  Real ***Diff_conf_fl = new Real**[3];
  for(int bz=0; bz<3; bz++){
    Diff_conf_fl[bz] = new Real*[Dim];
    for(int d=0; d<Dim; d++){
      Diff_conf_fl[bz][d] = new Real[Nt/2-1];
    }
  }
  Real ***Err_conf_fl = new Real**[3];
  for(int bz=0; bz<3; bz++){
    Err_conf_fl[bz] = new Real*[Dim];
    for(int d=0; d<Dim; d++){
      Err_conf_fl[bz][d] = new Real[Nt/2-1];
    }
  }
  
  for(int bz=0; bz<NBz; bz++){
    for(int d=0; d<Dim; d++){
      for(int t=0; t<Nt/2-1; t++){
	
	for(int m=0; m<NMass; m++){
	  Diff_conf_fl[bz][d][t] += e2/16*C[m]*C[m]*Media[bz][m][d][t];
	  Err_conf_fl[bz][d][t] += pow(C[m]*C[m]/16,2)*Err[bz][d][m][t];
	}
	Err_conf_fl[bz][d][t] = sqrt(Err_conf_fl[bz][d][t]);
      }
    }
  }
  
  
  Real ***Diff_conf_fl_d = new Real**[3];
  for(int bz=0; bz<3; bz++){
    Diff_conf_fl_d[bz] = new Real*[2];
    for(int d=0; d<2; d++){
      Diff_conf_fl_d[bz][d] = new Real[Nt/2-1];
    }
  }
  Real ***Err_conf_fl_d = new Real**[3];
  for(int bz=0; bz<3; bz++){
    Err_conf_fl_d[bz] = new Real*[2];
    for(int d=0; d<2; d++){
      Err_conf_fl_d[bz][d] = new Real[Nt/2-1];
    }
  }
  for(int bz=0; bz<NBz; bz++){
    for(int t=0; t<Nt/2-1; t++){
      
      for(int d=0; d<2; d++){
	Diff_conf_fl_d[bz][0][t] +=Diff_conf_fl[bz][d][t];
	Err_conf_fl_d[bz][0][t] += Err_conf_fl[bz][d][t]*Err_conf_fl[bz][d][t]/4;
      }
      Diff_conf_fl_d[bz][0][t] = Diff_conf_fl_d[bz][0][t]/2;
      Err_conf_fl_d[bz][0][t]=Err_conf_fl_d[bz][0][t]/4;
    }
  }

  FILE *Correlators_bz_outputs_M;
  char open_Correlators_bz_outputs_M[1024];

  for(int bz=0; bz<NBz; bz++){

    sprintf(open_Correlators_bz_outputs_M, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Differenze/T40L48_beta4.14_b%d/Diff_b%d_M_xy.out", Bz[bz], Bz[bz]);
    if ((Correlators_bz_outputs_M = fopen(open_Correlators_bz_outputs_M, "w")) == NULL ){
      printf("Error opening the input file: %s\n",open_Correlators_bz_outputs_M);
      exit(EXIT_FAILURE);
    }
    
    for(int t=0; t<Nt/2-1; t++){

      fprintf(Correlators_bz_outputs_M, "%d " "%s " "%s\n", t, conv(Diff_conf_fl_d[bz][0][t]).c_str(), conv(Err_conf_fl_d[bz][0][t]).c_str());
      
    }

  }
    
  /*
  Real ******Corr_bz_boot = new Real*****[3];
  for(int bz=0; bz<3; bz++){
    Corr_bz_boot[bz] = new Real****[NConf[1]];
    for(int iconf = 0; iconf < NConf[1]; iconf++){
      Corr_bz_boot[bz][iconf] = new Real***[NMass];
      for(int m = 0; m < NMass; m++){
	Corr_bz_boot[bz][iconf][m] = new Real**[Dim];
	for(int d=0; d<Dim; d++){
	  Corr_bz_boot[bz][iconf][m][d] = new Real*[Nt/2-1];
	  for(int t=0; t<Nt/2-1; t++){
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
	for(int t=0; t<Nt/2-1; t++){
	  for(int iboot=0; iboot<Nboot; iboot++){
	    for(int iconf=0; iconf<NConf[bz]; iconf++){
	      Corr_bz_boot[bz][iconf][m][d][t][iboot] = Diff[bz][distrib(gen)][m][d][t];
	    } 
	  }
	  
	}
      }
    }
  }
  
  
  
  
  for(int bz=0; bz<NBz; bz++){
    
    
    for(int iconf=0; iconf<NConf[bz]; iconf++){
      for(int d=0; d<Dim; d++){
	for(int m=0; m<NMass; m++){
	  for(int t=0; t<Nt/2-1; t++){
	    
	    cout << "bz: " << bz << " iconf: " << iconf << " d: " << d << " m: " << endl;
	    cout << "DIFF: t " << t << "  " <<  Boot_Mean(Corr_bz_boot[bz][iconf][m][d][t], Nboot) << "  " << Boot_Sigma(Corr_bz_boot[bz][iconf][m][d][t], Nboot) << endl;
	    
	  }
	}
      }
    }
    
  }
  */

  return 0;
  

}
