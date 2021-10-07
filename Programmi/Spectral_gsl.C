//g++ A-Wall -O0 -ggdb3 -lgmp -lgmpxx -std=c++11 -I/mnt/c/Users/navigl/Desktop/gsl/include -o Spectral_gsl Spectral_gsl.C `gsl-config --cflags --libs`

//g++ -o Spectral_gsl Spectral_gsl.C -O0 -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx -L/Users/manuel/Desktop/gmpfrxx  -g -Wall -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas

#include <fstream>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <sstream>
#include <vector>
#include <cmath>
#include "gmpfrxx.h"


using namespace std;
using T=mpfr_class;


/////////////////////////////////////////                                   
//SETTA LA PRECISIONE DESIDERATA IN BITS                                   
const int P = 560;
/////////////////////////////////////////


/////////////////////////////////////////
//GLOBAL QUANTITIES

int D_Latt = 10;
T Corr[10], t_in[10];
T t_glb_R;

//Parameters of the computation
T beta=1, bar_omega=0;
////////////////////////////////////////

string conv(const mpfr_class& in)
{
  ostringstream os;
  os.precision(mpf_get_default_prec()/4);
  
  os<<in;

  return os.str();
}

mpfr_class conv(const string& in)
{
  istringstream is(in);

  mpfr_class out;

  is>>out;

  return out;
}




struct params_t{
  
  T t_i;
  T t_j;
  T t;
  T beta;
  T bar_omega;
  params_t(T t_i, T t_j, T beta, T bar_omega) : t_i(t_i), t_j(t_j), beta(beta), bar_omega(bar_omega) {}
  
};


T K(T omega, T t, T beta){
  
  //cout << "KK om: " << omega << "  " << t << "  " << beta << endl;
  //cout << "KK: " << cosh(omega*(t - beta/2)) << "  " << sinh(beta*omega/2)<< endl;
  //cout << "KK: " << cosh(omega*(t - beta/2))/sinh(beta*omega/2)<< endl;
  
  T A=cosh(omega*(t - beta/2));
  T B=sinh(beta*omega/2);
  return A/B*omega; //Aggiunta f(\omega)=\omega necessaria per regolarizzazione in \omega=0
}

double Integrand_W(double s, void *arg){ 
  
  params_t *params=(params_t*)arg;
  
  T ti = params->t_i;
  T tj = params->t_j;
  T beta = params->beta;
  T bar_omega = params->bar_omega;
  
  
  
  T A;
  A =K(s,ti,beta)*(s-bar_omega)*(s-bar_omega)*K(s,tj,beta);
  double Int = A.get_d();
  
  
  return Int;
  
}


double Integrand_R(double s, void *arg){

  params_t *params=(params_t*)arg;

  T t = params->t;
  T beta = params->beta;
  t=t_glb_R;
  return K(s,t,beta).get_d();
  
}




int main(){
  
  //////////////// PASSO LA PRECISIONE SETTATA DI DEFAULT //////////////    
  mpfr_set_default_prec(P);



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
  for(int i=0; i<D_Latt; i++){
    
    fscanf(Correlators_Inputs, "%lf " "%s " "%lf\n", &trash1, trash3, &trash2);
    t_in[i] = trash1/10;
    Corr[i] = conv(trash3);
    cout << "t[" << i << "]=" << t_in[i] <<  " Corr[" << i << "]=" << Corr[i] << endl; 
    
  }
  

  fclose(Correlators_Inputs);

  //FINE INPUT DATA
  
  
  
  
  
  //PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA
  params_t params_Graph(0.1,0.2,beta,bar_omega);
  
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

  
  for(double i=0; i<1000; i++){
    
    fprintf(Integrand_beha, "%lf " "%s\n", i/1000+0.0100, conv(Integrand_W(i/1000+0.01, &params_Graph)).c_str());
    fprintf(K_beha, "%lf " "%s\n", i/1000+0.0100, conv(K(i/1000+0.0100,0.1, beta)).c_str());
    
  }

  fclose(Integrand_beha);
  fclose(K_beha);
  //FINE PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA




  
  //INIZIO INTEGRAZIONE
  
  int workspace_size=10000000;
  
  double t_i=0.1, t_j=0.1;
  
  
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  
  
  gsl_function f;
  params_t params(t_i,t_j,beta,bar_omega);
  f.function=Integrand_W;
  f.params=&params;
  
  
  double result, abserr;
  double start=0, epsabs=0,epsrel=1e-6;
  
  
  
  gsl_integration_qagiu(&f,start,epsabs,epsrel,workspace_size, workspace, &result, &abserr);

  cout << "W: " << result << endl;
  
  
  t_glb_R = t_i;
  f.function=Integrand_R;
  gsl_integration_qagiu(&f, start, epsabs,epsrel,workspace_size, workspace, &result, &abserr);
  
  
  gsl_integration_workspace_free(workspace);
  
  
  cout << "R: " << result << endl;;
  
  
  //FINE INTEGRAZIONE




  
  return 0;
  
}
