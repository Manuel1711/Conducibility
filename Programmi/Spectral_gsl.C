//g++ -std=c++11 -o Spectral_gsl Spectral_gsl.C -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include -L/usr/local/include  -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas


#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath>
#include "mp.h"
#include "smear.h"

/////////////////////////////////////////
//GLOBAL QUANTITIES
Real Corr[10], t_in[tmax];
Real t_glb_R;
/////////////////////////////////////////


int main(){
  
  cout << R_NInt(infLimit, supLimit, 3) << endl;
  cout << W_NInt(infLimit, supLimit, 3, 2) << endl;
  cout << f_NInt(infLimit, supLimit, 3) << endl;     
  
  //////////////// PASSO LA PRECISIONE SETTATA DI DEFAULT //////////////    
  PrecMatr W_Mat(tmax,tmax), Id(tmax, tmax), Id_bis(tmax, tmax);
  PrecVec R(tmax);
  
  
  
  //INPUT DATA

  //Input Correlatori
  FILE *Correlators_Inputs;
  char open_Correlators_Inputs[1024];
  
  sprintf(open_Correlators_Inputs, "Correlatori_Corrected/zero_mus/b3.85/corr_mu=0.000_format");
  
  if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs);
    exit(EXIT_FAILURE);
  }
  
  
  
  double trash1, trash2;
  char trash3[1024];
  for(int i=0; i<tmax; i++){
    
    fscanf(Correlators_Inputs, "%lf " "%s " "%lf\n", &trash1, trash3, &trash2);
    //t_in[i] = trash1;
    Corr[i] = conv(trash3);
    cout  <<  " Corr[" << i << "]=" << Corr[i] << endl; 
    
  }
  

  fclose(Correlators_Inputs);


  //Binnaggio t_in

  for(int i=1; i<tmax+1; i++){
    t_in[i-1] = i;
    //cout << "t_in[" << i-1 << "]" << t_in[i-1] << endl;
  }
    
  //FINE INPUT DATA
  
  
  
  
  //PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA
  params_t params_Graph(0.2,0.1,beta,Estar);

   
  
  FILE *Integrand_beha, *K_beha;
  char open_Integrand_beha[1024], open_K_beha[1024];
   
  sprintf(open_Integrand_beha, "Output/Integrand_beha.out");
  sprintf(open_K_beha, "Output/K_beha.out");
  
  if ((Integrand_beha = fopen(open_Integrand_beha, "w")) == NULL ){
    printf("Error opening the output file: %s\n",open_Integrand_beha);
    exit(EXIT_FAILURE);
  }
  
  if ((K_beha = fopen(open_K_beha, "w")) == NULL ){
    printf("Error opening the output file: %s\n",open_K_beha);
    exit(EXIT_FAILURE);
  }

  
  for(double i=0; i<100000; i++){
    
    //fprintf(Integrand_beha, "%lf " "%s\n", i/100, conv(Integrand_W(i/100, &params_Graph)).c_str());
    fprintf(K_beha, "%lf " "%s\n", i/100, conv(K(i/100,0.1, beta)).c_str());
    
  }

  fclose(Integrand_beha);
  fclose(K_beha);
  //FINE PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA
  


  
  // ************************ INIZIO METODO ******************************
  
  
  //CALCOLO MATRICE W E VETTORE R
  
  for(int i=0; i<tmax; i++){
    for(int j=0; j<tmax; j++){

#if defined(EXP)
      W_Mat(i,j) = W_an_exp(t_in[i], t_in[j], Estar);
#endif
#if defined(COS)
      W_Mat(i,j) = W_NInt(infLimit, supLimit, t_in[i], t_in[j]);
#endif
      cout << "W[" << t_in[i] << "][" << t_in[j] << "]=" << W_Mat(i,j) << endl;
      
    }//j
    
    
    //Calcolo vettore R
#if defined(HLN)
    R(i) = 1/(t_in[i])*exp(-E0*t_in[i]);
#endif
    
#if defined(BG)
#if defined(EXP)
    R(i) = 1/t_in[i];
#endif
#if defined(COS)
    R(i) = R_NInt(infLimit, supLimit, t_in[i]);
#endif
#endif
    
    cout << "R[" << t_in[i] << "]=" << R(i) << endl;
	 
    
  }//i
  
  
  

  
  
  // INVERSIONE MATRICE W
  
  const auto Winv=W_Mat.inverse();
  
  //Id = ((Wm1*W_Mat)-Eigen::Identity(31,31)).norm();
  for(int i=0; i<tmax; i++){
    for(int j=0; j<tmax; j++){
      
      //cout << "Id: " << Id(i,j) << endl;
      
    }
  } 
  // FINE INVERSIONE MATRICE W
  
  
  //CALCOLO f
  PrecVec f(tmax);
#if defined(HLN)
  for(int i=0; i<tmax; i++){
#if defined(EXP)
    f(i) = N(t_in[i])*D(t_in[i]);
#endif
#if defined(COS)
    f(i) = f_NInt(infLimit, supLimit,t_in[i]);
#endif
    cout << "f: " << f(i) << endl;
  }
#endif
  // FINE CALCOLO f
 

  // CALCOLO g
  Real den =  R.transpose()*Winv*R;
  PrecVec g;
#if defined(HLN)
  Real numA = R.transpose()*Winv*f;
  Real num = 1-numA;
  PrecVec g1 = Winv*f;
  g=Winv*f+ Winv*R*num/den;
#endif
  
#if defined(BG)
  g=Winv*R/den;
#endif
  // FINE CALCOLO g


  //Output coefficienti
  for(int i=0; i<tmax; i++) cout << "g: " << g(i) << endl; 
    
  FILE *q_t_out;
  char open_q_t_out[1024];
  
  sprintf(open_q_t_out, "Output/q_t_out.out");
  
  if ((q_t_out = fopen(open_q_t_out, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_q_t_out);
    exit(EXIT_FAILURE);
  }
  
  
  for(int i=0; i<tmax; i++){
    
    fprintf(q_t_out, "%s " "%s\n", conv(t_in[i]).c_str(), conv(g(i)).c_str());
    
  }
  
  
  fclose(q_t_out);
  
  

  //Output funzione di smearing
  FILE *Delta_S;
  char open_Delta_S[1024];
  
  sprintf(open_Delta_S, "Output/Delta_Smear.out");
  
  if ((Delta_S = fopen(open_Delta_S, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Delta_S);
    exit(EXIT_FAILURE);
  }
  
  
#if defined(BG)
    E0=0;
#endif
    fprintf(Delta_S, "@type xy\n");
    for(double i=0; i<300; i++){

      fprintf(Delta_S, "%s " "%s\n", conv(E0 +i/100).c_str(), conv(Delta_Smear(E0 + i/100, g, t_in)).c_str());
      
    }
    
#if defined(HLN)
    fprintf(Delta_S, "\n \n @type xy \n");
    
    for(double i=0; i<300; i++){
      
      fprintf(Delta_S, "%s " "%s\n", conv(E0 +i/100).c_str(), conv(Target_F(E0 + i/100)).c_str());
    }
  
    
    fprintf(Delta_S, "\n \n @type xy \n");
    for(double i=0; i<300; i++){
      Real df=Target_F(E0 + i/100)-Delta_Smear(E0 + i/100, g, t_in);
      fprintf(Delta_S, "%s " "%s\n", conv(E0 +i/100).c_str(), conv(df).c_str());
    }
    
    
    fclose(Delta_S);
    
    
    //Output differenza target/HLN
    FILE *Diff;
    char open_Diff[1024];
    
    
    sprintf(open_Diff, "Output/Diff.out");
    
    if ((Diff = fopen(open_Diff, "w")) == NULL ){
      printf("Error opening the input file: %s\n",open_Diff);
      exit(EXIT_FAILURE);
    }
    
    fprintf(Diff, "\n \n @type xy \n");
    for(double i=0; i<300; i++){
      Real df=Target_F(E0 + i/100)-Delta_Smear(E0 + i/100, g, t_in);
      fprintf(Diff, "%s " "%s\n", conv(E0 +i/100).c_str(), conv(df).c_str());
    }
    
    fclose(Diff);
#endif

    
    return 0;
   
}
