#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <cmath>
#include "mp.h"
#include "smear.h"

int main(){
  
  PrecVec Corr(Nt), q(Nt);
  Real t_in[Nt];
  

  
  //INPUT QUANTITIES


  //Input Correlatori
  FILE *Correlators_Inputs;
  char open_Correlators_Inputs[1024];
  
  //sprintf(open_Correlators_Inputs, "Correlatori_Corrected/zero_mus/b3.85/corr_mu=0.000_format");
  sprintf(open_Correlators_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/Jack_Corr/Corr_Jack_41_xx");
  
  if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs);
    exit(EXIT_FAILURE);
  }
  double trash1, trash2;
  char trash3[1024];
  for(int i=0; i<Nt; i++){
    
    //fscanf(Correlators_Inputs, "%lf " "%s " "%lf\n", &trash1, trash3, &trash2);
    fscanf(Correlators_Inputs, "%lf " "%s\n" , &trash1, trash3);
    Corr(i) = conv(trash3);
    cout  <<  " Corr[" << i << "]=" << Corr(i) << endl; 
    
  }
  
  
  fclose(Correlators_Inputs);



  //Input Coefficients
  FILE *Coefficients_Inputs;
  char open_Coefficients_Inputs[1024];
  
  //sprintf(open_Correlators_Inputs, "Correlatori_Corrected/zero_mus/b3.85/corr_mu=0.000_format");
  sprintf(open_Coefficients_Inputs, "/Users/manuel/Documents/GitHub/Conducibility/Programmi/Output/q_t_out.out");
  
  if ((Coefficients_Inputs = fopen(open_Coefficients_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Coefficients_Inputs);
    exit(EXIT_FAILURE);
  }
  
  char trash4[1024], trash5[1024];
  for(int i=0; i<Nt; i++){
    fscanf(Coefficients_Inputs, "%s " "%s\n" , trash5, trash4);
    t_in[i] = conv(trash5);
    q(i) = conv(trash4);
    cout << "q: " << q(i) << endl;
  }
  
  
  //Prova
  PrecVec Corr_try(Nt);
  for(int i=0; i<Nt; i++){
    const auto f_try=
      [=](const Real& E) -> Real
      { 
#if defined(EXP)
       return exp(-E)*K(E, t_in[i], beta);
#endif
#if defined(COS)
       return exp(-E)*K(E, t_in[i], beta)/E;  
#endif
      };
    
    Corr_try(i)=
      bq::gauss_kronrod<Real,61>::integrate(f_try,infLimit,supLimit,5,1e-16);
    //cout << "Try: " << Corr_try(i) << endl;
  }
  
  
  //FINE INPUT QUANTITIES
  
  
  
  
  // Spectral function computation
#if defined(EXP)
  Real fomega =1;
#endif
#if defined(COS)
  Real fomega = Estar;
#endif
  
  cout << "rho(" << Estar << ")=" << Estar*spectral(q, Corr_try) << endl;
  cout << "rho_true(" << Estar << ")=" << exp(-Estar) << endl;
#if defined(HLN)
  cout << "rho_int(" << Estar << ")=" << rho_NInt(infLimit, supLimit, Estar, sigma) << endl;
#endif
  
  
  
  
  
  
  
  
  return 0;
  
  
}
