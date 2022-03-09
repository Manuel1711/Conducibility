#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath> 
#include "mp.h"
#include "smear.h"
#include "statistical.h"


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
  PrecMatr W_Mat(Nt,Nt), Id(Nt, Nt), Id_bis(Nt, Nt), Corr_e(Nboot,Nt), Corr_o(Nboot,Nt);
  PrecVec Corr_err_e(Nt), Corr_err_o(Nt), Corr_mu_e(Nt), Corr_mu_o(Nt);
  Real rho[Nboot], rho_S[Nboot], MIN=0;
  
  
  
  // INPUT QUANTITIES
  //Input Correlatori
  FILE *Correlators_Inputs, *Correlators_Inputs_MuS;
  char open_Correlators_Inputs[1024], open_Correlators_Inputs_MuS[1024];
  
  sprintf(open_Correlators_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Bootstraps/T40L48_beta4.14_b93_Subtracted/boot_correlators_xy.out");
  sprintf(open_Correlators_Inputs_MuS, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Means_Sigmas/T40L48_beta4.14_b93_Subtracted/mu_sigma_correlators_xy.out");
  //sprintf(open_Correlators_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Bootstraps/T20L24_beta3.787_b93_Subtracted/boot_correlators_z.out");
  //2sprintf(open_Correlators_Inputs_MuS, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Our_Correlators_bz/Subtracted/Means_Sigmas/T20L24_beta3.787_b93_Subtracted/mu_sigma_correlators_z.out");
  
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
  int e=0,o=0;
  int em=0, om=0;
  for(int i=0; i<2*Nt+1; i++){ //C'è +1 perché invece di t=0 prendiamo t=20 (parte da 0 perché deve scorrere).
    for(int iboot=0; iboot<Nboot; iboot++){
      fscanf(Correlators_Inputs, "%lf " "%s\n", &trash1, trash2);
      if(i%2 != 0){
	Corr_o(iboot, o) = conv(trash2);
        cout  <<  " Corr_o[" << o << "," << iboot << "]=" << Corr_o(iboot,o) << endl;
	if(iboot==Nboot-1) o++;
      }
      else if(i%2==0 and i>0){ 
      	Corr_e(iboot,e) = conv(trash2);
	cout  <<  " Corr_e[" << e << "," << iboot << "]=" << Corr_e(iboot,e) << endl; 
	if(iboot==Nboot-1) e++;
      }
     
    }//iboot

    double trash;
    char trash_b[1024], trash_bb[1024];
    fscanf(Correlators_Inputs_MuS, "%lf " "%s " "%s\n", &trash, trash_b, trash_bb);
    if(i%2 != 0 and i>0){
      Corr_mu_o(om) = conv(trash_b);
      Corr_err_o(om) = conv(trash_bb);
      om++;
      cout << "ERR: " << Corr_err_o(em) << " i: " << i << endl;
    }
    else if(i%2==0 and i>0){ 
      Corr_mu_e(em) = conv(trash_b);
      Corr_err_e(em) = conv(trash_bb);
      cout << "ERR: " << Corr_err_e(em) << " i: " << i << endl;
      em++;
    }
    
    
  }//i
  e=0;
  o=0;
  em=0;
  om=0;
  
  fclose(Correlators_Inputs);



  
  
  
  // ************************ INIZIO METODO ******************************
  
  FILE *Lambda_Shape_out;
  char open_Lambda_Shape_out[1024];

  sprintf(open_Lambda_Shape_out, "Output/Lambda_Shape.out");
  if ((Lambda_Shape_out = fopen(open_Lambda_Shape_out, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Lambda_Shape_out);
    exit(EXIT_FAILURE);
  }

  
  for(int ilambda=40; ilambda<100; ilambda++){

    Real lambda = conv(to_string(ilambda))/100;
    cout << "Lambda: " << lambda << endl;
     
    for(int iboot=0; iboot<Nboot; iboot++){
      //Inputs
      PrecVec Corr(Nt), Corr_err(Nt), Corr_Mu(Nt), R(Nt);
      PrecMatr Cov(Nt,Nt);
      Real t_a[Nt];
      for(int i=0; i<Nt; i++){
	if(EO==0){
	  Corr(i)=Corr_e(iboot,i);
	  Corr_err(i) = Corr_err_e(i);
	  Corr_Mu(i) = Corr_mu_e(i);
	}
	if(EO==1){
	  Corr(i)=Corr_o(iboot,i);
	  Corr_err(i) = Corr_err_o(i);
	  Corr_Mu(i) = Corr_mu_o(i);
	}
	//cout << "iboot: " << iboot << " i: " << i << " C: " << Corr(i) << endl;
      } 
      
      //Binnaggio t
      e=0;
      o=0;
      for(int i=1; i<Nt+1; i++){
	if(EO==0){
	  t_a[e] = 2*i;
	  e++;
	}
	if(EO==1){
	  t_a[o] = 2*i-1; 
	  o++;
	}
	//cout << "t_a[" << i-1 << "]" << t_a[i-1] << endl;
      }
      for(int i=0; i<Nt; i++){
	for(int j=0; j<Nt; j++){
	  if(i==j){
	    Cov(i,j)=Corr_err(i)*Corr_err(j);
	    //cout << "Cov: " <<  Cov(i,j) << endl; 
	  }
	  else Cov(i,j)=0;
	}
      }
      
      
      
      //Calcolo matrice W
      for(int i=0; i<Nt; i++){
	for(int j=0; j<Nt; j++){
#if defined(EXP)
	  W_Mat(i,j) = (1-lambda)*W_an_exp(t_a[i], t_a[j], Estar)+lambda*Cov(i,j);
#endif 
#if defined(COS)
	  W_Mat(i,j) = (1-lambda)*W_NInt(infLimit, supLimit, t_a[i], t_a[j], Estar)+lambda*Cov(i,j);
	  //cout << "i=" << i << " j=" << j << " A: " << W_NInt(infLimit, supLimit, t_a[i], t_a[j], Estar) << " B: " << Cov(i,j)  << endl;
	  //cout << "LAMBDA:  i" << i << " j=" << j << " A: " << (1-lambda)*W_NInt(infLimit, supLimit, t_a[i], t_a[j], Estar) << " B: " << lambda*Cov(i,j)  << endl;
#endif
	  //cout << "W[" << t_a[i] << "][" << t_a[j] << "]=" << W_Mat(i,j) << endl;
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
	//cout << "R[" << t_a[i] << "]=" << R(i) << endl;
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
      for(int i=0; i<Nt; i++)
	//cout << "g(" << i << ")="<< g(i) << endl;
      
      
      //Output coefficienti
      /*FILE *q_t_out;
	char open_q_t_out[1024];
	
	sprintf(open_q_t_out, "Output/q_t_out.out");
	if ((q_t_out = fopen(open_q_t_out, "w")) == NULL ){
	printf("Error opening the output file: %s\n",open_q_t_out);
	exit(EXIT_FAILURE);
	}
	for(int i=0; i<Nt; i++){
	fprintf(q_t_out, "%s " "%s\n", conv(t_a[i]).c_str(), conv(g(i)).c_str());
	}
	fclose(q_t_out);
      */
      
      
      //Output funzione di smearing
      /*FILE *Delta_S;
      char open_Delta_S[1024];
      
      if(EO==0)sprintf(open_Delta_S, "Output/Delta_Smear_e.out");
      if(EO==1)sprintf(open_Delta_S, "Output/Delta_Smear_o.out");
      
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
      */

      
      // Spectral function computation
#if defined(EXP)
      Real fomega =1;
#endif
#if defined(COS)
      Real fomega = Estar;
#endif
      rho[iboot]=spectral(g, Corr);
      rho_S[iboot] = stat_unc(g, Corr_err);
      
      if(abs(rho[iboot])>MIN) MIN=abs(rho[iboot]);
      
      //cout << "Bg ****** rho(" << iboot << ")=" << rho(iboot,eo) << endl;
      if(EO==0) cout << "rho_even(" << Estar << ")=" << rho[iboot] << "   " << stat_unc(g, Corr_err) << "  Lambda=" << lambda << endl;
      if(EO==1) cout << "rho_odd(" << Estar << ")=" << rho[iboot] << "   " << stat_unc(g, Corr_err) << " Lambda=" << lambda << endl;
      //cout << "rho_true(" << Estar << ")=" << exp(-Estar) << endl;
#if defined(HLN)
      //cout << "rho_int(" << Estar << ")=" << rho_NInt(infLimit, supLimit, Estar, sigma) << endl;
#endif
      
      
    }//iboot
    
    //cout << "MIN=" << MIN << endl;
    
    //cout << "Sigma: " << rho(0)+rho(1) << endl;
    //cout << "Sigma Plot: " << (rho(0)+rho(1))/(2*(4*Pi/(137.04)*(0.6666666*0.6666666 + 2*0.33333333*0.33333333))) << endl;
    
    
    //Media e Sigma bootstrap
    Real rho_mu, rho_sigma, Bar_Sigma; 
    
    for(int iboot=0; iboot<Nboot; iboot++){
      //cout << "rho(" << iboot << ")=" << rho[iboot] << " pm " << rho_S[iboot] << endl;
    }
    
    rho_mu = Boot_Mean(rho, Nboot);
    rho_sigma = Boot_Sigma(rho, Nboot);
    
    if(EO==0) cout << "rho_even_MU(" << Estar << ")=" << rho_mu << "   " << rho_sigma << endl;
    if(EO==1) cout << "rho_odd_MU(" << Estar << ")=" << rho_mu << "   " << rho_sigma << endl;
    
    fprintf(Lambda_Shape_out, "%s " "%s " "%s\n", conv(lambda).c_str(), conv(rho_mu).c_str(), conv(rho_sigma).c_str());
    
    
    
  }//lambda
  
  fclose(Lambda_Shape_out);
  
  /*for(int iboot=0; iboot<Nboot; iboot++){
    rho_mu += rho(iboot);
    cout << "rho(" << iboot << ")=" << rho(iboot) << endl;
    Bar_Sigma += rho(iboot)*rho(iboot);
    //cout << "rho: " << rho(iboot,eo) << endl;
  }
  Bar_Sigma = Bar_Sigma/Nboot;
  rho_mu = rho_mu/Nboot;
  //cout << "AAA: " << Bar_Sigma(eo) << "  " << rho_mu(eo)*rho_mu(eo) << endl;
  rho_sigma = sqrt((Bar_Sigma-rho_mu*rho_mu)/Nboot);*/
  
  
  //Real Sigma_mu = rho_mu(0)+rho_mu(1);
  //Real Sigma_s = sqrt(rho_sigma(0)*rho_sigma(0) + rho_sigma(1)*rho_sigma(1));
  
  //cout << "Sigma: " << Sigma_mu  << "   " << Sigma_s << endl;
  //cout << "Sigma Plot: " << beta*Sigma_mu/(2*(4*Pi/(137.04)*(0.6666666*0.6666666 + 2*0.33333333*0.33333333))) <<  "  " <<  Sigma_s/(2*(4*Pi/(137.04)*(0.6666666*0.6666666 + 2*0.33333333*0.33333333))) << endl;
  
  
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
