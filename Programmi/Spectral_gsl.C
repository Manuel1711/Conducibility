//g++ -Wall -O0 -ggdb3 -lgmp -lgmpxx -std=c++11 -I/mnt/c/Users/navigl/Desktop/gsl/include -o Spectral_gsl Spectral_gsl.C `gsl-config --cflags --libs`

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


int K_choice = 1;

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
T beta=10, bar_omega=5;
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
params_t(T t, T beta, T bar_omega) : t(t), beta(beta), bar_omega(bar_omega) {}
  
  
};


T K(T omega, T t, T beta){
  
  //cout << "KK om: " << omega << "  " << t << "  " << beta << endl;
  //cout << "KK: " << cosh(omega*(t - beta/2)) << "  " << sinh(beta*omega/2)<< endl;
  //cout << "KK: " << cosh(omega*(t - beta/2))/sinh(beta*omega/2)<< endl;
  
  T ret;
  
  if(K_choice == 1){
    
    ret= exp(-omega*t) + exp(-(beta-t)*omega);
    
  }
  
  else if(K_choice == 0){
    
    T A=cosh(omega*(t - beta/2));
    T B=sinh(beta*omega/2);
    ret = A/B*omega;
    
  }

  
  return ret; //Aggiunta f(\omega)=\omega necessaria per regolarizzazione in \omega=0
}

double Integrand_W(double s, void *arg){ 
  
  params_t *params=(params_t*)arg;
  
  T ti = params->t_i;
  T tj = params->t_j;
  T beta = params->beta;
  T bar_omega = params->bar_omega;

  

  
  
  T A;
  
  A = K(s,ti,beta)*(s-bar_omega)*(s-bar_omega)*K(s,tj,beta);
  //cout << "A: " << A << endl;
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



T q_i(T **W, T R[], int i){
  
  T N=0, D=0;
  for(int j=1; j<D_Latt; j++){
    
    N += 1/W[i][j]*R[j];
    
    
    
    for(int k=1; k<D_Latt; k++){
      
      D += R[k]*1/W[k][j]*R[j];
      
    }
  }
  
  return N/D;

}



T Delta_Smear(T omega, T q[], T t_in[]){
  
  T D;
  for(int i=0; i<D_Latt; i++){
    
    //cout << "q[i]: " << q[i] << " K: " << K(omega, t_in[i], beta) << "  " << t_in[i] <<  endl;
    
    D += q[i]*K(omega, t_in[i], beta);
    
  }
  
  return D;
  
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
    t_in[i] = trash1;
    Corr[i] = conv(trash3);
    cout << "t[" << i << "]=" << t_in[i] <<  " Corr[" << i << "]=" << Corr[i] << endl; 
    
  }
  

  fclose(Correlators_Inputs);
  
  //FINE INPUT DATA
  
  
  
  
  //PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA
  params_t params_Graph(0.2,0.1,beta,bar_omega);

  
  
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
    
    fprintf(Integrand_beha, "%lf " "%s\n", i/100, conv(Integrand_W(i/100, &params_Graph)).c_str());
    fprintf(K_beha, "%lf " "%s\n", i/100, conv(K(i/100,0.1, beta)).c_str());
    
  }

  fclose(Integrand_beha);
  fclose(K_beha);
  //FINE PLOT INTEGRANDI DI W E R IN FUNZIONE DI OMEGA
  
  
  
  
  //INIZIO INTEGRAZIONE
  
  int workspace_size=10000000;
  T** W = new T*[D_Latt];
  for(int i=0; i<D_Latt; i++){
    
    W[i] = new T[D_Latt];
    
  }
  
  T R[D_Latt];
  
  
  
  //T t_i=0.1, t_j=0.1;

  
  /*FILE *q_beha;
  char open_q_beha[1024];
  sprintf(open_q_beha, "Output/q_beha.out");

  if ((q_beha = fopen(open_q_beha, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_q_beha);
    exit(EXIT_FAILURE);
  }
  
  for(int k=0; k<50; k++){
  
  */
  
  
  for(int i=1; i<D_Latt; i++){
    
    
    
      
    gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
    
    
    
    gsl_function f; 
    
    
    double abserr;
    double start=0, epsabs=0,epsrel=1e-10;
    
    for(int j=1; j<D_Latt; j++){
      
      
      double res_temp;
      params_t params(t_in[i], t_in[j], beta,bar_omega);
      f.function=Integrand_W;
      f.params=&params;
      
      
      
      gsl_integration_qagiu(&f,start,epsabs,epsrel,workspace_size, workspace, &res_temp, &abserr);
      
      W[i][j] = res_temp;
      cout << "W[" << i << "][" << j << "]=" << W[i][j] << endl;
      
      
    }//j
    
    
    
    double res_temp2;
    
    t_glb_R = t_in[i];
    params_t params_R(t_glb_R, beta,bar_omega);
    f.function=Integrand_R;
    f.params=&params_R;
    gsl_integration_qagiu(&f, start, epsabs,epsrel,workspace_size, workspace, &res_temp2, &abserr);
    
    R[i] = res_temp2;
    
    gsl_integration_workspace_free(workspace);
    
    //cout << "RES: " << res_temp2;
    cout << "R[" << i << "]=" << R[i] << endl;;
    
    
    
  }//i
  
  
  
  //FINE INTEGRAZIONE
  
  
  
  T q[D_Latt];
  
  
  for(int i=1; i<D_Latt; i++){
    
    
    q[i] = q_i(W, R, i);
    cout << q_i(W, R, i) << endl;;
    
  }
  
  /*
    fprintf(q_beha, "%d " "%lf\n", k, q_i(W, R, 2).get_d());
    
    
    }//k
    
    
    fclose(q_beha);
  */
  
  
  FILE *q_t_out;
  char open_q_t_out[1024];

  sprintf(open_q_t_out, "Output/q_t_out.out");
  
  if ((q_t_out = fopen(open_q_t_out, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_q_t_out);
    exit(EXIT_FAILURE);
  }
  

  for(int i=0; i<D_Latt; i++){
    
    fprintf(q_t_out, "%d " "%lf\n", i, q[i].get_d());
    
  }
  
  
  fclose(q_t_out);
  
  
  
  FILE *Delta_S;
  char open_Delta_S[1024];
  
  sprintf(open_Delta_S, "Output/Delta_Smear.out");
  
  if ((Delta_S = fopen(open_Delta_S, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Delta_S);
    exit(EXIT_FAILURE);
  }
  
  
  

  
  for(double i=1; i<400; i++){
    
    if(i<200) fprintf(Delta_S, "%lf " "%lf\n", (-200+i)/10, Delta_Smear((-200+i)/10, q, t_in).get_d());
    else fprintf(Delta_S, "%lf " "%lf\n", (i-200)/10, Delta_Smear((i-200)/10, q, t_in).get_d());
    
  }
  
  fclose(Delta_S);
  
  return 0;
  
}
