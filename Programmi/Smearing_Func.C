#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath> 
#include "mp.h"
#include "smear.h"

/////////////////////////////////////////
//GLOBAL QUANTITIES
Real t_a[Nt];
/////////////////////////////////////////


int main(){
  
#if defined(SINGLE)
  Real Estar=0.0;
#endif
#if defined(CICLO)
  FILE *Estar_Inputs;
  char open_Estar_Inputs[1024];
  
  sprintf(open_Estar_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/bin.out");
  
  if ((Estar_Inputs = fopen(open_Estar_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Estar_Inputs);
    exit(EXIT_FAILURE);
  }
  char tt[1024];
  fscanf(Estar_Inputs, "%s", tt);
  
  Real Estar = conv(tt);
  
#endif
 
  cout << "************ Estar=" << Estar << " ************" << endl;
  
  //////////////// PASSO LA PRECISIONE SETTATA DI DEFAULT //////////////    
  PrecMatr W_Mat(Nt,Nt), Id(Nt, Nt), Id_bis(Nt, Nt);
  PrecVec R(Nt), Corr(Nt), Corr_e(Nt), Corr_o(Nt), Corr_err(Nt), Corr_err_e(Nt), Corr_err_o(Nt), rho(2);
  
  
  
  // INPUT QUANTITIES
  
  //Input Correlatori
  FILE *Correlators_Inputs;
  char open_Correlators_Inputs[1024];
  
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Desktop/Russi/1/b3.85/corr_mu=0.140_format");
  sprintf(open_Correlators_Inputs, "/Users/manuel/Desktop/zero_must/b3.85/corr_mu=0.285_format");
  //sprintf(open_Correlators_Inputs, "Correlatori_Corrected/zero_mus/b3.85/corr_mu=0.000_format");
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/Jack_Corr/Corr_Jack_41_xx");
  
  if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs);
    exit(EXIT_FAILURE);
  }
  double trash1;
  char trash3[1024], trash2[1024];
  for(int i=0; i<2*Nt; i++){
    
    fscanf(Correlators_Inputs, "%lf " "%s " "%s\n", &trash1, trash3, trash2);
    //fscanf(Correlators_Inputs, "%lf " "%s\n" , &trash1, trash3);
    if(i%2){}
    else{
      cout << "i: " << i << endl;
      Corr_e(i/2) = conv(trash3);
      Corr_err_e(i/2) = conv(trash2);
      cout  <<  " Corr_e[" << i/2 << "]=" << Corr_e(i/2) << endl; 
    }
    if(i%2){
      cout << "i: " << i << endl;
      Corr_o(i/2) = conv(trash3);
      Corr_err_o(i/2) = conv(trash2);
      cout  <<  " Corr_o[" << i/2 << "]=" << Corr_o(i/2) << endl;
    }
  }
  
  
  fclose(Correlators_Inputs);
  
  //FOR SULLE PARI E LE DISPARI
  for(int eo=0; eo<2; eo++){






    
    if(eo==0){
      Corr=Corr_e;
      Corr_err = Corr_err_e;
    }
    if(eo==1){
      Corr=Corr_o;
      Corr_err = Corr_err_o;
    }
      //Binnaggio t
  for(int i=0; i<Nt; i++){
#if defined(EV)
    t_a[i-tmin] = 2*i;
#endif
#if defined(ODD)
    t_a[i-tmin] = 2*i+1; 
#endif
    cout << "t_a[" << i-tmin << "]" << t_a[i-tmin] << endl;
  }
  

  // FINE INPUT QUANTITIES


  
  
  //Prova
  /*
  PrecVec Corr_try(Nt);
  for(int i=0; i<Nt; i++){
    const auto f_try=
      [=](const Real& E) -> Real
      { 
#if defined(EXP)
       return exp(-E)*K(E, t_a[i], beta);
#endif
#if defined(COS)
       return exp(-E)*K(E, t_a[i], beta)/E;  
       #endif
      };
    
    Corr_try(i)=
      bq::gauss_kronrod<Real,61>::integrate(f_try,infLimit,supLimit,5,1e-16);
    //cout << "Try: " << Corr_try(i) << endl;
  }
  */

  
  //PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA
  
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
    fprintf(K_beha, "%lf " "%s\n", i/100, conv(K(i/100,0.00988*2, 0.00988*10)).c_str());
    
  }
  
  fclose(Integrand_beha);
  fclose(K_beha);
  //FINE PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA
  
  
  
  
  // ************************ INIZIO METODO ******************************
  
  
  //CALCOLO MATRICE W E VETTORE R
  
  for(int i=0; i<Nt; i++){
    for(int j=0; j<Nt; j++){
#if defined(EXP)
      W_Mat(i,j) = W_an_exp(t_a[i], t_a[j], Estar);
#endif
#if defined(COS)
      W_Mat(i,j) = W_NInt(infLimit, supLimit, t_a[i], t_a[j], Estar);
#endif
      cout << "W[" << t_a[i] << "][" << t_a[j] << "]=" << W_Mat(i,j) << endl;
      
    }//j
      
    
    //Calcolo vettore R
#if defined(HLN)
#if defined(EXP)
    R(i) = 1/(t_a[i])*exp(-E0*t_a[i]);
#endif
#if defined(COS)
    R(i) = R_NInt(infLimit, supLimit, t_a[i]);
#endif
#endif
    
#if defined(BG)
#if defined(EXP)
    R(i) = 1/(t_a[i]);
#endif
#if defined(COS)
    R(i) = R_NInt(infLimit, supLimit, t_a[i]);
#endif
#endif
    
    cout << "R[" << t_a[i] << "]=" << R(i) << endl;
    
  }//i
  
  
  
  
  // INVERSIONE MATRICE W
  const auto Winv=W_Mat.inverse();
  
  
  //CALCOLO f
  PrecVec f(Nt);
#if defined(HLN)
  f = f_func(t_a, sigma, Estar);
#endif
  
  
  // CALCOLO g
  PrecVec g;
  g = Coeff(R,Winv,f);
  
  

  
  //Output coefficienti
  for(int i=0; i<Nt; i++) cout << "g: " << g(i) << endl; 
  
  FILE *q_t_out;
  char open_q_t_out[1024];
  
  sprintf(open_q_t_out, "Output/q_t_out.out");
  
  if ((q_t_out = fopen(open_q_t_out, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_q_t_out);
    exit(EXIT_FAILURE);
  }
  
  
  for(int i=0; i<Nt; i++){
    
    fprintf(q_t_out, "%s " "%s\n", conv(t_a[i]).c_str(), conv(g(i)).c_str());
    
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
    
    fprintf(Delta_S, "%s " "%s\n", conv(E0 +0.0001 + i/100).c_str(), conv(Delta_Smear(E0 + 0.0001 + i/100, g, t_a)).c_str());
      
  }
  
#if defined(HLN)
  fprintf(Delta_S, "\n \n @type xy \n");
  
  for(double i=0; i<300; i++){
    
    fprintf(Delta_S, "%s " "%s\n", conv(E0 + 0.0001 +i/100).c_str(), conv(Target_F(E0 + 0.0001 + i/100, Estar, sigma)).c_str());
  }
  
  
  fprintf(Delta_S, "\n \n @type xy \n");
  for(double i=0; i<300; i++){
    Real df=Target_F(E0 + 0.0001 + i/100, Estar, sigma)-Delta_Smear(E0 + 0.0001 + i/100, g, t_a);
    fprintf(Delta_S, "%s " "%s\n", conv(E0 + 0.0001 + i/100).c_str(), conv(df).c_str());
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
    Real df=Target_F(E0 + 0.0001 + i/100, Estar, sigma)-Delta_Smear(E0 + 0.0001 + i/100, g, t_a);
    fprintf(Diff, "%s " "%s\n", conv(E0 + 0.0001 + i/100).c_str(), conv(df).c_str());
  }
  
  fclose(Diff);
#endif
  
  // Spectral function computation
#if defined(EXP)
  Real fomega =1;
#endif
#if defined(COS)
  Real fomega = Estar;
#endif
  
  if(eo==0) cout << "rho_even(" << Estar << ")=" << spectral(g, Corr) << "   " << stat_unc(g, Corr_err) << endl;
  if(eo==1) cout << "rho_odd(" << Estar << ")=" << spectral(g, Corr) << "   " << stat_unc(g, Corr_err) << endl;
  //cout << "rho_true(" << Estar << ")=" << exp(-Estar) << endl;
#if defined(HLN)
  //cout << "rho_int(" << Estar << ")=" << rho_NInt(infLimit, supLimit, Estar, sigma) << endl;
#endif
  


  }//eo



  
#if defined(CICLO)
  
  FILE *rho_output;
  char open_rho_output[1024];

  sprintf(open_rho_output, "/tmp/temp.out");
  if ((rho_output = fopen(open_rho_output, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_rho_output);
    exit(EXIT_FAILURE);
  }

  fprintf(rho_output,"%s " "%s\n", conv(spectral(g, Corr)).c_str(), conv(abs(stat_unc(g, Corr_err))).c_str());
  
  fclose(rho_output);
  
#endif

  
  // Procedura delta_sigma
  
  //delta_sigma_procedure(0.05, Estar, t_a, R, Winv);
  
  

  
  return 0;
    
} 
