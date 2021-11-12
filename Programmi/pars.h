//////////////////////////////
//PARAMETRI DEL CALCOLO
static T beta=10, Estar=0.5;
#define D_Latt 15
//Nazario shifta di 1. Per lui t_max=30 partendo in realtà da 0 (Quindi D_Latt=31). Qui si parte sempre da 1.
static const T sigma=0.1;
T E0=0.1; 
static const T alpha=0;
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

int K_choice = 1;


/////////////////////////////////////////                                   
//SETTA LIMITE INFERIORE E SUPERIORE INTEGRAZIONE NUMERICA
const Real infLimit=0.23532;
const Real supLimit=inf;			
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
