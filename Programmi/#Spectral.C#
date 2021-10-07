#include <fstream>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
//#include <TH1F.h>

using namespace std;

double Corr[10];


double K(double omega, double t, double beta){
  
  return cosh(omega*(t - beta/2))/sinh(beta*omega/2);
  
}


double Integrand(double s, double ti, double tj, double beta, double bar_omega){ //Trasformazione argomento integrale per rendere l'intervallo finito
  
  //cout << "ti: " << ti <<  " tj: " << tj << " s: " << s  << " s/(1-s) " << s/(1-s) <<  " beta: " << beta << " K(xi, s/(1-s)): " << K(s/(1-s), ti, beta) << endl;


  //L'integrazione in d\omega diventa in ds
  
  return K(s/(1-s), ti, beta)*pow((s/(1-s)-bar_omega),2)*K(s/(1-s), tj, beta)/(1-s*s);
  
  
}


double Integral_TM(int Nbins, double x_min, double x_max, double ti, double tj, double beta, double bar_omega){
  
  double step = (x_max-x_min)/Nbins;
  double x=x_min;
  double inte=0;
  cout << "KKK " << inte << endl;
  
  while(x<0.9){
    
    //cout << "KKK " << inte << endl;
    //cout << "KK " << inte << endl;
    //cout << " KKK " << Integrand(x, ti, tj, beta, bar_omega) + Integrand(x+step, ti, tj, beta, bar_omega) << endl;
    
    double sum = Integrand(x, ti, tj, beta, bar_omega) + Integrand(x+step, ti, tj, beta, bar_omega);
    
    inte = inte + sum*0.5*step;
    
    cout << "SUM: " << inte << endl;
    
    
    x += step;
    
    
    
  }
  
  return inte;
  
}



int main(){
  

  //Input Correlatori
  FILE *Correlators_Inputs;
  char open_Correlators_Inputs[1024];
  
  
  sprintf(open_Correlators_Inputs, "Correlators/zero_mus/b3.85/corr_mu=0.000_format");

  if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs);
    exit(EXIT_FAILURE);
  }
  
  
  double trash1, trash2;
  for(int i=0; i<10; i++){
    
    fscanf(Correlators_Inputs, "%lf " "%lf " "%lf\n", &trash1, &Corr[i], &trash2);
    cout << "Corr[" << i << "] " << Corr[i] << endl; 
    
  }

  cout << K(2,3,4) << endl;


  double W_try; // provo a calcolare una entry di W (senza correlazione) per bar_omega = 0 e lambda=1

  W_try = Integral_TM(10000, 0.01, 1, 1, 2, 0.002, 0); //Ho messo x_min=0.01 perchè per \omega=0 K ha una singolarità
  
  
  cout << "W: " << W_try << endl;
  
  return 0;

}
