#ifndef _MP_H
#define _MP_H
#endif

Real Boot_Mean(Real Boot[], int Nboot){
  Real Mean=0;
  for(int iboot=0; iboot<Nboot; iboot++)
    Mean += Boot[iboot];
  return Mean/Nboot;
}

/*Real Boot_Sigma(Real Boot[], int Nboot, Real Mean){
  
  Real Mean_Sq=0;
  for(int iboot=0; iboot<Nboot; iboot++)
    Mean_Sq += Boot[iboot]*Boot[iboot];
  Mean_Sq = Mean_Sq/Nboot;

  return sqrt(Mean_Sq-Mean*Mean);
  
}
*/
Real Boot_Sigma(Real Boot[], int Nboot){
  
  Real Mean_Sq=0, Mean=0;
  for(int iboot=0; iboot<Nboot; iboot++){
    Mean_Sq += Boot[iboot]*Boot[iboot];
    Mean += Boot[iboot];
  }
  Mean_Sq = Mean_Sq/Nboot;
  Mean = Mean/Nboot;
  
  return sqrt(Mean_Sq-Mean*Mean);
  
}
