#include "pars.h"

 
T Z(){

  return (1+erf(Estar/(sqrt(2)*sigma)))/2;
  
}


T Target_F(T E){
  
   return exp(-pow(E-Estar,2)/(2*sigma*sigma))/(sqrt(2*Pi)*sigma*Z());
  
}

T W_an(T ti, T tj, T Estar){
  
  //return -(-2+2*Estar*(ti+tj)-pow(Estar,2)*pow(ti+tj,2))/(pow(ti+tj,3));
#if defined(HLN)
  return exp(-(ti+tj-alpha)*E0)/(ti+tj-alpha);
#endif
#if defined(BG)
  return -(-2+2*Estar*(ti+tj)-pow(Estar,2)*pow(ti+tj,2))/pow(ti+tj,3);
#endif
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
  
#if defined(EXP)
    
    ret= exp(-omega*t); //+ exp(-(beta-t)*omega);
    
#endif
  
#if defined(COS)
    
    T A=cosh(omega*(t - beta/2));
    T B=sinh(beta*omega/2);
    ret = A/B;
    
#endif
  
  return ret;
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
