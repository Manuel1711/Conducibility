#include <fstream>
#include <iostream>
#include <math.h>
//#include <TH1F.h>

using namespace std;

double Corr[10];


double b_T(double omega, double t, double beta){

  return cosh(beta - t/2)/sinh(beta*omega/2);
  
}


int main(){
  

  //Input Correlatori
  FILE *Correlators_Inputs;
  char open_Correlators_Inputs[1024];
  
  
  sprintf(open_Correlators_Inputs, "Correlators/zero_mus/b3.85/corr_mu=0.000_format");

  if ((Correlators_Inputs = fopen(open_Correlators_Inputs, "r")) == NULL ){
    printf("Error opening the input file: %s\n",open_Correlators_Inputs);
    exit(EXIT_FAILURE);
  }
  
  
  double trash1, trash2;
  for(int i=0; i<10; i++){
    
    fscanf(Correlators_Inputs, "%lf " "%lf " "%lf\n", &trash1, &Corr[i], &trash2);
    cout << "Corr[" << i << "] " << Corr[i] << endl; 
    
  }

  cout << b_T(2,3,4) << endl;
  
  return 0;

}
