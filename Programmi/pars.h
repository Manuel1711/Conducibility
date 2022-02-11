//////////////////////////////
//PARAMETRI DEL CALCOLO
#define tmax 39
#define tmin 0
#define Nt 10 
static Real beta=tmax+1;
//static Real beta=1/200;
//Nazario shifta di 1. Per lui t_max=30 partendo in realtà da 0 (Quindi D_Latt=31). Qui si parte sempre da 1.
static const Real sigma=1.7/12;
static Real E0=0.0; 
static const Real alpha=0;
static const Real lambda=0.0;
static const int Nboot = 100;
/////////////////////////////////////////                                   
 
//SETTA LA PRECISIONE DESIDERATA IN BITS                                   
const int P = 1024;
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
#if defined(BG)
const Real infLimit=0;
#endif
#if defined(HLN)
const Real infLimit=E0;
#endif
const Real supLimit=inf;			
/////////////////////////////////////////

struct params_t{
  
  Real t_i;
  Real t_j;
  Real t;
  Real beta;
  Real bar_omega;
params_t(Real t_i, Real t_j, Real beta, Real bar_omega) : t_i(t_i), t_j(t_j), beta(beta), bar_omega(bar_omega) {}
params_t(Real t, Real beta, Real bar_omega) : t(t), beta(beta), bar_omega(bar_omega) {}
  
  
};
