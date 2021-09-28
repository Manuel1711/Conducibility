//g++ -Wall -I/mnt/c/Users/navigl/Desktop/gsl/include -o Spectral Spectral.C `gsl-config --cflags --libs`

#include <fstream>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>


using namespace std;

double Corr[10];

struct params_t{
  
  double t_i;
  double t_j;
  double t;
  double beta;
  double bar_omega;
  params_t(double t_i, double t_j, double beta, double bar_omega) : t_i(t_i), t_j(t_j), beta(beta), bar_omega(bar_omega) {}
  
};


double K(double omega, double t, double beta){
  
  
  return cosh(omega*(t - beta/2))/sinh(beta*omega/2);
  
}
 
double Integrand(double s, void *arg){ 
  
  params_t *params=(params_t*)arg;
  
  double ti = params->t;
  double tj = params->t;
  double beta = params->beta;
  double bar_omega = params->bar_omega;
  
  return K(s,ti,beta)*pow((s-bar_omega),2)*K(s,tj,beta);
  
}








int main(){

  double t_i=1, t_j=1, beta=0.01, bar_omega=0;
  params_t params(t_i,t_j,beta,bar_omega);
  
  //Plot integrando in funzione di s
  FILE *Integrand_beha;
  char open_Integrand_beha[1024];

  sprintf(open_Integrand_beha, "Output/Integrand_beha.out");
  
  if ((Integrand_beha = fopen(open_Integrand_beha, "w")) == NULL ){
    printf("Error opening the output file: %s\n",open_Integrand_beha);
    exit(EXIT_FAILURE);
  }

  for(double i=0; i<100; i++){

    fprintf(Integrand_beha, "%lf " "%lf\n", i/10, Integrand(i, &params));

  }

  fclose(Integrand_beha);
  


  //Integrale
  int workspace_size=100000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  
  
  gsl_function f;
  f.function=Integrand;
  f.params=&params;
  
  
  double result, abserr;
  double start=0.2, epsabs=0,epsrel=1e-6;
  
  gsl_integration_qagiu(&f,start,epsabs,epsrel,workspace_size, workspace, &result, &abserr);
  
  
  gsl_integration_workspace_free(workspace);
  
  return 0;
  
}
