#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath> 
#include "mp.h"
#include "smear.h"



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
  PrecMatr W_Mat(Nt,Nt), Id(Nt, Nt), Id_bis(Nt, Nt), Corr_e(Nboot,Nt), Corr_o(Nboot,Nt), rho(Nboot,2);
  PrecVec R(Nt), Corr(Nt), Corr_err(Nt), Corr_err_e(Nt), Corr_err_o(Nt);
  
  
  
  // INPUT QUANTITIES
  //Input Correlatori
  FILE *Correlators_Inputs, *Correlators_Inputs_MuS;
  char open_Correlators_Inputs[1024], open_Correlators_Inputs_MuS[1024];
  
  //sprintf(open_Correlators_Inputs_mIz, "/Users/manuel/Desktop/Russi/1/b3.85/corr_mu=0.000_format");
  //sprintf(open_Correlators_Inputs_mIz, "/Users/manuel/Desktop/zero_must/b3.93998000/corr_33_mu=0.0");
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Desktop/zero_must/b3.93998000/corr_mu=0.140_format");
  sprintf(open_Correlators_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Bootstraps/T20L24_beta3.787_b93_Subtracted/boot_correlators_xy.out");
  sprintf(open_Correlators_Inputs_MuS, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Means_Sigmas/T20L24_beta3.787_b93_Subtracted/mu_sigma_correlators_xy.out");
  //sprintf(open_Correlators_Inputs, "Correlatori_Corrected/zero_mus/b3.85/corr_mu=0.000_format");
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/Jack_Corr/Corr_Jack_41_xx");
  if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs);
    exit(EXIT_FAILURE);
  }
  if ((Correlators_Inputs_MuS = fopen(open_Correlators_Inputs_MuS, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs_MuS);
    exit(EXIT_FAILURE);
  }
  double trash1;
  char trash2[1024];
  for(int i=0; i<2*Nt; i++){
    for(int iboot=0; iboot<Nboot; iboot++){
      
      fscanf(Correlators_Inputs, "%lf " "%s\n", &trash1, trash2);
      if(i%2){}
      else{
	Corr_e(iboot,i/2) = conv(trash2);
	cout  <<  " Corr_e[" << i/2 << "]=" << Corr_e(iboot,i/2) << endl; 
      }
      if(i%2){
	Corr_o(iboot, i/2) = conv(trash2);
        cout  <<  " Corr_o[" << i/2 << "]=" << Corr_o(iboot,i/2) << endl;
      }
    }//iboot

    if(i%2){} 
    else{
      double trash;
      char trash_b[1024], trash_bb[1024];
      fscanf(Correlators_Inputs_MuS, "%lf " "%s " "%s\n", &trash, trash_b, trash_bb);
      Corr_err_e(i/2) = conv(trash_bb);
    }
    if(i%2){
      double trash;
      char trash_b[1024], trash_bb[1024];
      fscanf(Correlators_Inputs_MuS, "%lf " "%s " "%s\n", &trash, trash_b, trash_bb);
      Corr_err_o(i/2) = conv(trash_bb);
    }
    
  }//i
  
  
  fclose(Correlators_Inputs);



  
  
  
  // ************************ INIZIO METODO ******************************
  
  
  //Inputs
  for(int eo=0; eo<2; eo++){//pari e dispari
    for(int iboot=0; iboot<Nboot; iboot++){
      
      PrecMatr Cov(Nt,Nt);
      Real t_a[Nt];
      cout << " OOOOOOO " << endl;
      for(int i=0; i<Nt; i++){
	if(eo==0){
	  Corr(i)=Corr_e(iboot,i);
	  Corr_err(i) = Corr_err_e(i);
	}
	if(eo==1){
	  Corr(i)=Corr_o(iboot,i);
	  Corr_err(i) = Corr_err_o(i);
	}
      }
      cout << " OOOOOOO222 " << endl;
      //Binnaggio t
      for(int i=0; i<Nt; i++){
	if(eo==0)
	  t_a[i-tmin] = 2*i;
	if(eo==1)
	  t_a[i-tmin] = 2*i+1; 
	cout << "t_a[" << i-tmin << "]" << t_a[i-tmin] << endl;
      }
      for(int i=0; i<Nt; i++){
	for(int j=0; j<Nt; j++){
	  if(i==j){
	    Cov(i,j)=Corr_err(i)*Corr_err(j);
	  }
	  else Cov(i,j)=0;
	}
      }
  
  
  
  //Calcolo matrice W
  for(int i=0; i<Nt; i++){
    for(int j=0; j<Nt; j++){
#if defined(EXP)
      W_Mat(i,j) = (1-lambda)*W_an_exp(t_a[i], t_a[j], Estar)+(lambda)/pow(Corr(0),2)*Cov(i,j);
#endif 
#if defined(COS)
      W_Mat(i,j) = (1-lambda)*W_NInt(infLimit, supLimit, t_a[i], t_a[j], Estar)+(lambda)/pow(Corr(0),2)*Cov(i,j);
      cout << "PPP: " << Corr(0) << endl;
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
  
  
  
  
  //Inversione matrice W
  const auto Winv=W_Mat.inverse();
  
  
  //Calcolo f (solo BG modificato)
  PrecVec f(Nt);
#if defined(HLN)
  f = f_func(t_a, sigma, Estar);
#endif
  
  
  //Calcolo g
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
  
  if(eo==0)sprintf(open_Delta_S, "Output/Delta_Smear_e.out");
  if(eo==1)sprintf(open_Delta_S, "Output/Delta_Smear_o.out");
  
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
  
  rho(iboot,eo)=spectral(g, Corr); 
  if(eo==0) cout << "rho_even(" << Estar << ")=" << rho(iboot,eo) << "   " << stat_unc(g, Corr_err) << endl;
  if(eo==1) cout << "rho_odd(" << Estar << ")=" << rho(iboot,eo) << "   " << stat_unc(g, Corr_err) << endl;
  //cout << "rho_true(" << Estar << ")=" << exp(-Estar) << endl;
#if defined(HLN)
  //cout << "rho_int(" << Estar << ")=" << rho_NInt(infLimit, supLimit, Estar, sigma) << endl;
#endif
  

    }//iboot
  }//eo

  

  //cout << "Sigma: " << rho(0)+rho(1) << endl;
  //cout << "Sigma Plot: " << (rho(0)+rho(1))/(2*(4*Pi/(137.04)*(0.6666666*0.6666666 + 2*0.33333333*0.33333333))) << endl;

  
  //Media e Sigma bootstrap
  PrecVec rho_mu(2), rho_sigma(2); 
  
  for(int eo=0; eo<2; eo++){
    PrecVec Bar_Sigma(2);
    for(int iboot=0; iboot<Nboot; iboot++){
      rho_mu(eo) += rho(iboot,eo);
      Bar_Sigma(eo) += rho(iboot,eo)*rho(iboot,eo);
    }
    rho_mu(eo) = rho_mu(eo)/Nboot;
    rho_sigma(eo) = (Bar_Sigma(eo)-rho_mu(eo)*rho_mu(eo)/Nboot)/Nboot;
  } 
  
  cout << "Sigma: " << rho_mu(0)+rho_mu(1) << "   " << sqrt(rho_sigma(0)*rho_sigma(0) + rho_sigma(1)*rho_sigma(1)) << endl;
  cout << "Sigma Plot: " << (rho_mu(0)+rho_mu(1))/(2*(4*Pi/(137.04)*(0.6666666*0.6666666 + 2*0.33333333*0.33333333))) << endl;
  
  
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
