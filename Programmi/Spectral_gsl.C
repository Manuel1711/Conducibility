//g++ -Wall -O0 -ggdb3 -lgmp -lgmpxx -std=c++11 -I/mnt/c/Users/navigl/Desktop/gsl/include -o Spectral_gsl Spectral_gsl.C `gsl-config --cflags --libs`

//g++ -o Spectral_gsl Spectral_gsl.C -O0 -I /usr/include/eigen3/ -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx -L/Users/manuel/Desktop/gmpfrxx  -g -Wall -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas

//g++ -o Spectral_gsl Spectral_gsl.C -O0 -I /usr/include/eigen3/ -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include/boost/ -L/usr/local/lib  -g -Wall -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas

//g++ -std=c++11 -o Spectral_gsl Spectral_gsl.C -I/mnt/c/Users/navigl/Desktop/gsl/include -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include -L/usr/local/include  -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas


#include <fstream>
#include <iostream>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <sstream>
#include <vector>
#include <cmath>
#include "gmpfrxx.cpp" 
#include <eigen3/Eigen/Dense>
#include "mp.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/multiprecision/mpfr.hpp>

namespace bm=boost::multiprecision;
namespace bq=boost::math::quadrature;
 
using Real=
  bm::number<bm::mpfr_float_backend<128>>;

/////////////////////////////////////////                                   
//SETTA LA PRECISIONE DESIDERATA IN BITS                                   
const int P = 512;

struct Initer
{
  Initer()
  {
    mpfr_class::set_dprec(P);
  }
};
 
Initer initer;
/////////////////////////////////////////


/////////////////////////////////////////                                   
//SETTA LIMITE INFERIORE E SUPERIORE INTEGRAZIONE NUMERICA
const Real inf=
  std::numeric_limits<Real>::infinity();
const Real infLimit=0.23532;
const Real supLimit=inf;			
/////////////////////////////////////////


#define Pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

using namespace std;
using namespace Eigen;
using T=mpfr_class;

using PrecMatr= Matrix<T,Dynamic,Dynamic>;
using PrecVec = Matrix<T,Dynamic, 1>;
 

int K_choice = 1;

//Scelta metodo
// 0=BG, 1=Nazario
int Method_choice = 1;


//Parameters of the computation
T beta=10, Estar=0.5;
#define D_Latt 16
//Nazario shifta di 1. Per lui t_max=30 partendo in realtà da 0 (Quindi D_Latt=31). Qui si parte sempre da 1.
T sigma=0.1;
T E0=0.1;
T alpha=0;
////////////////////////////////////////


/////////////////////////////////////////
//GLOBAL QUANTITIES

T Corr[10], t_in[D_Latt];
T t_glb_R;


/////////////////////////////////////////

struct params_t{
  
  T t_i;
  T t_j;
  T t;
  T beta;
  T bar_omega;
params_t(T t_i, T t_j, T beta, T bar_omega) : t_i(t_i), t_j(t_j), beta(beta), bar_omega(bar_omega) {}
params_t(T t, T beta, T bar_omega) : t(t), beta(beta), bar_omega(bar_omega) {}
  
  
};



//INIZIO FUNZIONI


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



T Z(){

  return (1+erf(Estar/(sqrt(2)*sigma)))/2;
  
}


T Target_F(T E){
  
   return exp(-pow(E-Estar,2)/(2*sigma*sigma))/(sqrt(2*Pi)*sigma*Z());
  
}

T W_an(T ti, T tj, T Estar){
  
  //return -(-2+2*Estar*(ti+tj)-pow(Estar,2)*pow(ti+tj,2))/(pow(ti+tj,3));
  if(Method_choice==1) return exp(-(ti+tj-alpha)*E0)/(ti+tj-alpha);
  else if(Method_choice==0) return -(-2+2*Estar*(ti+tj)-pow(Estar,2)*pow(ti+tj,2))/pow(ti+tj,3);
}




T N(T t){
  
  //cout << " N: " << 1/(2*Z())*exp(((alpha-t)*((alpha-t)*pow(sigma,2)+2*Estar))/2) << endl;
  return 1/(2*Z())*exp(((alpha-t)*((alpha-t)*pow(sigma,2)+2*Estar))/2);
}

T D(T t){
  
  cout << "D: " << 1+erf(((alpha-t)*pow(sigma,2)+Estar-E0)/(sqrt(2)*sigma)) << " erf: " << erf(((alpha-t)*pow(sigma,2)+Estar-E0)/(sqrt(2)*sigma)) <<"   " << (alpha-t)*pow(sigma,2)+Estar-E0   << endl;
  return  1+erf(((alpha-t)*pow(sigma,2)+Estar-E0)/(sqrt(2)*sigma)); 
  
}


T K(T omega, T t, T beta){
  
  
  T ret;
  
  if(K_choice == 1){
    
    ret= exp(-omega*t); //+ exp(-(beta-t)*omega);
    
  }
  
  else if(K_choice == 0){
    
    T A=cosh(omega*(t - beta/2));
    T B=sinh(beta*omega/2);
    ret = A/B*omega;
    
  }

  
  return ret; //Aggiunta f(\omega)=\omega necessaria per regolarizzazione in \omega=0
}




T q_i(T **W, T R[], int i){
  
  T N=0, D=0;
  for(int j=0; j<D_Latt; j++){
    
    N += W[i][j]*R[j];
    
    
    
    for(int k=0; k<D_Latt; k++){
      
      D += R[k]*W[k][j]*R[j];
      
    }
  }
  
  return N/D;

}



T Delta_Smear(T omega, PrecVec q, T t_in[]){
  
  T D;
  for(int i=0; i<D_Latt; i++){
    
    //cout << "q[i]: " << q(i) << " K: " << K(omega, t_in[i], beta) << " t: " << t_in[i] << " omega: " << omega <<  endl;
    
    D += q(i)*K(omega, t_in[i], beta);
    
  }
  
  return D;
  
}

// FINE FUNZIONI





int main(){

  

  
  //////////////// PASSO LA PRECISIONE SETTATA DI DEFAULT //////////////    
  PrecMatr W_Mat(D_Latt,D_Latt), Id(D_Latt, D_Latt), Id_bis(D_Latt, D_Latt);
  PrecVec R(D_Latt);
  
  
  
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
    //t_in[i] = trash1;
    Corr[i] = conv(trash3);
    cout  <<  " Corr[" << i << "]=" << Corr[i] << endl; 
    
  }
  

  fclose(Correlators_Inputs);


  //Binnaggio t_in

  for(int i=1; i<D_Latt+1; i++){
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
  
  
  //CALCOLO MATRICE W

  for(int i=0; i<D_Latt; i++){
    
    
    for(int j=0; j<D_Latt; j++){
      
      W_Mat(i,j) = W_an(t_in[i], t_in[j], Estar);
      //cout << "Analitico ---> W[" << t_in[i] << "][" << t_in[j] << "]=" << W_Mat(i,j) << endl;
      
    }//j
    
    
    if(Method_choice==1)R(i) = 1/(t_in[i])*exp(-E0*t_in[i]);
    else if(Method_choice==0)R(i) = 1/t_in[i];
      
    cout << "R[" << t_in[i] << "]=" << R(i) << endl;
	 
    
  }//i
  
  //FINE CALCOLO MATRICE W
  

  
  
  // INVERSIONE MATRICE W
  
  
  cout << "W_matrix: " << endl;
  FILE *W_out;
  char open_W_out[1024];
  
  sprintf(open_W_out, "Output/WM_out.out");
  
  if ((W_out = fopen(open_W_out, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_W_out);
    exit(EXIT_FAILURE);
  }
  
  for(int col=0; col<D_Latt; col++){
    for(int row=0; row < D_Latt ; row++){
      
      fprintf(W_out, "%lf ", W_Mat(col,row).get_d());
      cout << W_Mat(col,row) << "  ";
      
    }

    fprintf(W_out, "\n");
    cout << endl;
    
  } 
  
  
  
  fclose(W_out);
  

  
  const auto Winv=W_Mat.inverse();
  
  /*cout << "W_matrix_inverse: " << endl;
  for(int col=0; col<D_Latt; col++){
    for(int row=0; row < D_Latt ; row++){
      
      cout << Winv(col,row) << "  ";
      
    }
    
    cout << endl;
    
    }*/
  
  //Id = ((Wm1*W_Mat)-Eigen::Identity(31,31)).norm();
  for(int i=0; i<D_Latt; i++){
    for(int j=0; j<D_Latt; j++){
      
      //cout << "Id: " << Id(i,j) << endl;
      
    }
  } 
  
  // FINE INVERSIONE MATRICE W


  PrecVec f(D_Latt);
  if(Method_choice==1){
  //CALCOLO f

  for(int i=0; i<D_Latt; i++){
    
     f(i) = N(t_in[i])*D(t_in[i]);
    cout << "f: " << f(i) << "  " << N(t_in[i]) << "  " << D(t_in[i]) <<  endl;
  }
  
  // FINE CALCOLO f
  }//if

  
  // CALCOLO g
  T den =  R.transpose()*Winv*R;
  PrecVec g;
  if(Method_choice==1){
    T numA = R.transpose()*Winv*f;
    T num = 1-numA;
    PrecVec g1 = Winv*f;
    g=Winv*f+ Winv*R*num/den;
  }
  else if(Method_choice==0){
    g=Winv*R/den;
  }
  // FINE CALCOLO g


  for(int i=0; i<D_Latt; i++) cout << "g: " << g(i) << endl; 
    
  FILE *q_t_out;
  char open_q_t_out[1024];
  
  sprintf(open_q_t_out, "Output/q_t_out.out");
  
  if ((q_t_out = fopen(open_q_t_out, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_q_t_out);
    exit(EXIT_FAILURE);
  }
  
  
  for(int i=0; i<D_Latt; i++){
    
    fprintf(q_t_out, "%lf " "%lf\n", t_in[i].get_d(), g(i).get_d());
    
  }
  
  
  fclose(q_t_out);
  
  
  
  FILE *Delta_S;
  char open_Delta_S[1024];
  
  sprintf(open_Delta_S, "Output/Delta_Smear.out");
  
  if ((Delta_S = fopen(open_Delta_S, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Delta_S);
    exit(EXIT_FAILURE);
  }
  
  
  
  if(Method_choice==0) E0=0;
  fprintf(Delta_S, "@type xy\n");
  for(double i=0; i<300; i++){

    fprintf(Delta_S, "%lf " "%lf\n", E0.get_d() +i/100, Delta_Smear(E0 + i/100, g, t_in).get_d());
    
  }

  if(Method_choice==1){
  fprintf(Delta_S, "\n \n @type xy \n");
  
  for(double i=0; i<300; i++){
    
    fprintf(Delta_S, "%lf " "%lf\n", E0.get_d() +i/100, Target_F(E0 + i/100).get_d());
  }

 
  fprintf(Delta_S, "\n \n @type xy \n");
  for(double i=0; i<300; i++){
    T df=Target_F(E0 + i/100)-Delta_Smear(E0 + i/100, g, t_in);
    fprintf(Delta_S, "%lf " "%lf\n", E0.get_d() +i/100, df.get_d());
  }
  
  
  fclose(Delta_S);


  FILE *Diff;
  char open_Diff[1024];


  sprintf(open_Diff, "Output/Diff.out");

  if ((Diff = fopen(open_Diff, "w")) == NULL ){
    printf("Error opening the input file: %s\n",open_Diff);
    exit(EXIT_FAILURE);
  }

  fprintf(Diff, "\n \n @type xy \n");
  for(double i=0; i<300; i++){
    T df=Target_F(E0 + i/100)-Delta_Smear(E0 + i/100, g, t_in);
    fprintf(Diff, "%lf " "%lf\n", E0.get_d() +i/100, df.get_d());
  }

  fclose(Diff);
  }//if
  
  return 0;
  
}
