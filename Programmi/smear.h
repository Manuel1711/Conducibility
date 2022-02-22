#include "pars.h"

Real Z( Real s, Real Es){
  return (1+erf(Es/(sqrt(2)*s)))/2;
}


Real Target_F(Real E, Real Es, Real s){
  //return exp(-pow(E-Es,2)/(2*s*s))/(sqrt(2*Pi)*s*Z(s,Es));
  return 1/pow((s*Pi/2),2)*E/(sinh(E/s)); 
}


Real W_an_exp(Real ti, Real tj, Real Es){
#if defined(HLN)
  return exp(-(ti+tj-alpha)*E0)/(ti+tj-alpha);
#endif
#if defined(BG)
  return -(-2+2*Es*(ti+tj)-pow(Es,2)*pow(ti+tj,2))/pow(ti+tj,3);
#endif
} 




Real N(Real t, Real s, Real Es){
  return 1/(2*Z(s,Es))*exp(((alpha-t)*((alpha-t)*pow(s,2)+2*Es))/2);
}

Real D(Real t, Real s, Real Es){
  return  1+erf(((alpha-t)*pow(s,2)+Es-E0)/(sqrt(2)*s)); 
}

Real K(Real omega, Real t, Real beta){
  
  Real ret;
  
#if defined(EXP)
  ret=exp(-omega*t); //+ exp(-(beta-t)*omega);
#endif
  
#if defined(COS)
  Real A=cosh(omega*(t - beta/2));
  Real B=sinh(beta*omega/2);
  ret = A/B*omega; //Se moltiplico per omega alla fine ottengo rho/omega che è quello che mi interessa
#endif 
  
  return ret;
}




Real f_NInt(Real infLimit, Real supLimit, Real ti, Real s, Real Es){
  
  const auto f_cs=
    [=](const Real& E) -> Real
    { 
     return K(E, ti, beta)*Target_F(E,Es,s);
    };
  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
  
}




PrecVec f_func(Real t_a[], Real s, Real Es){
  
  PrecVec f(Nt);
  
  for(int i=0; i<Nt; i++){
#if defined(EXP)
    cout << "t: " << t_a[i] << endl;
    f(i) = N(t_a[i], s, Es)*D(t_a[i], s, Es);
#endif
#if defined(COS)
    f(i) = f_NInt(infLimit, supLimit,t_a[i], s, Es);
#endif
    //cout << "f: " << f(i) << endl;
  }
  
  return f;
  
}

 


PrecVec Coeff(PrecVec R, PrecMatr Winv, PrecVec f){
  
  Real den =  R.transpose()*Winv*R;
  //cout << "den: " << den << endl;
  PrecVec g;
  
  
#if defined(HLN)
  Real numA = R.transpose()*Winv*f;
  Real num = 1-numA;
  PrecVec g1 = Winv*f;
  return Winv*f+ Winv*R*num/den;
#endif
  
  
#if defined(BG)
  return Winv*R/den;
#endif
  
}





Real Delta_Smear(Real omega, PrecVec q, Real t_in[]){
  
  Real D;
  for(int i=0; i<Nt; i++)
    D += q(i)*K(omega, t_in[i], beta);
  
  return D;
  
}




Real R_NInt(Real infLimit, Real supLimit, Real ti){
  
  const auto f_cs=
    [=](const Real& E) -> Real
    {
     return K(E, ti, beta);
    };
  
  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
  
}


Real W_NInt(Real infLimit, Real supLimit, Real ti, Real tj, Real Es){
  
  const auto f_cs=
    [=](const Real& E) -> Real
    {
#if defined(BG)
     return K(E, ti, beta)*(E - Es)*(E - Es)*K(E,tj,beta);
#endif
#if defined(HLN)
     return K(E, ti, beta)*K(E,tj,beta);
#endif
    };
  
  
  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
}





Real rho_NInt(Real infLimit, Real supLimit, Real Es, Real s){
  
  const auto f_cs=
    [=](const Real& E) -> Real
    { 
     return Es*exp(-E)*Target_F(E,Es,s)/E;
    };
  
  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
  
}



Real spectral(PrecVec q, PrecVec C){
  
  Real rho=0;
  for(int i=0; i<Nt; i++){
    rho += q(i)*C(i);
    cout << "Bg: " << q(i) << "  " << C(i) << "  " << q(i)*C(i) << "  " << rho << endl;
  }
  
  return 2*Pi*rho;
  
}


Real stat_unc(PrecVec q, PrecVec dC){

  Real err=0;
  for(int i=0; i<Nt; i++)
    err += q(i)*dC(i);
  return err;
  
}
 

void delta_sigma_procedure(Real delta_sigma, Real Es, Real t_a[], PrecVec R, PrecMatr Winv){
  
  Real LHS = abs(delta_sigma-1);
  Real sigma_f = 0.00001;
  Real step = 0.0001;
  cout << "LHS: " << LHS << endl;
  
  while(abs(LHS-Delta_Smear(Es, Coeff(R,Winv,f_func(t_a, sigma_f, Es)), t_a)/Target_F(Es, Es, sigma_f)) > 10*step){
    sigma_f += step;
    cout << "Diff: " << LHS - Delta_Smear(Es, Coeff(R,Winv,f_func(t_a, sigma_f, Es)), t_a)/Target_F(Es, Es, sigma_f) << "  Sigma: " << sigma_f << endl;
  }
  cout << "Sigma Finale: " << sigma_f << endl;
  
}



