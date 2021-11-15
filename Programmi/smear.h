#include "pars.h"


Real Z(){

  return (1+erf(Estar/(sqrt(2)*sigma)))/2;
  
}


Real Target_F(Real E){
  
   return exp(-pow(E-Estar,2)/(2*sigma*sigma))/(sqrt(2*Pi)*sigma*Z());
  
}

Real W_an_exp(Real ti, Real tj, Real Estar){
  
  //return -(-2+2*Estar*(ti+tj)-pow(Estar,2)*pow(ti+tj,2))/(pow(ti+tj,3));
#if defined(HLN)
  return exp(-(ti+tj-alpha)*E0)/(ti+tj-alpha);
#endif
#if defined(BG)
  return -(-2+2*Estar*(ti+tj)-pow(Estar,2)*pow(ti+tj,2))/pow(ti+tj,3);
#endif
}




Real N(Real t){
  
  //cout << " N: " << 1/(2*Z())*exp(((alpha-t)*((alpha-t)*pow(sigma,2)+2*Estar))/2) << endl;
  return 1/(2*Z())*exp(((alpha-t)*((alpha-t)*pow(sigma,2)+2*Estar))/2);
}

Real D(Real t){
  
  cout << "D: " << 1+erf(((alpha-t)*pow(sigma,2)+Estar-E0)/(sqrt(2)*sigma)) << " erf: " << erf(((alpha-t)*pow(sigma,2)+Estar-E0)/(sqrt(2)*sigma)) <<"   " << (alpha-t)*pow(sigma,2)+Estar-E0   << endl;
  return  1+erf(((alpha-t)*pow(sigma,2)+Estar-E0)/(sqrt(2)*sigma)); 
  
}


Real K(Real omega, Real t, Real beta){
  
  
  Real ret;
  
#if defined(EXP)
    
    ret= exp(-omega*t); //+ exp(-(beta-t)*omega);
    
#endif
  
#if defined(COS)
    
    Real A=cosh(omega*(t - beta/2));
    Real B=sinh(beta*omega/2);
    ret = A/B;
    
#endif
  
  return ret;
}

  

Real q_i(Real **W, Real R[], int i){
  
  Real N=0, D=0;
  for(int j=0; j<tmax; j++){
    
    N += W[i][j]*R[j];
    
    
    
    for(int k=0; k<tmax; k++){
      
      D += R[k]*W[k][j]*R[j];
      
    }
  }
  
  return N/D;

}



Real Delta_Smear(Real omega, PrecVec q, Real t_in[]){
  
  Real D;
  for(int i=0; i<tmax; i++){
    
    //cout << "q[i]: " << q(i) << " K: " << K(omega, t_in[i], beta) << " t: " << t_in[i] << " omega: " << omega <<  endl;
    
    D += q(i)*K(omega, t_in[i], beta);
    //cout << "D: " << q(i) << "  " << K(omega, t_in[i], beta) << "  " << exp(-omega*t_in[i]) << "  " <<  q(i)*K(omega, t_in[i], beta) << endl;
    
  }

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


Real W_NInt(Real infLimit, Real supLimit, Real ti, Real tj){

  const auto f_cs=
    [=](const Real& E) -> Real
    {
#if defined(BG)
     return K(E, ti, beta)*(E - Estar)*(E - Estar)*K(E,tj,beta);
#endif
#if defined(HLN)
     return K(E, ti, beta)*K(E,tj,beta);
#endif
    };
  

  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
}



Real f_NInt(Real infLimit, Real supLimit, Real ti){
  
  const auto f_cs=
    [=](const Real& E) -> Real
    { 
     return K(E, ti, beta)*Target_F(E);
    };

  return
    bq::gauss_kronrod<Real,61>::integrate(f_cs,infLimit,supLimit,5,1e-16);
  
}
