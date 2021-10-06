// compile with: g++ -o test test.cpp -lgsl -lcblas

#include <stdio.h>
#include <gsl/gsl_integration.h>

struct params_t
{
  double A;
  double omega;
  double x0;
  params_t(double A,double omega,double x0) : A(A),omega(omega),x0(x0) {}
};

//! kernel to be integrated numerically
double kernel(double x,void *arg)
{
  params_t *params=(params_t*)arg;
  double A=params->A;
  double omega=params->omega;
  double x0=params->x0;
  
  return A*sin(omega*(x-x0));
}

//! used to check
double kernel_primitive(double x,void *arg)
{
  params_t *params=(params_t*)arg;
  double A=params->A;
  double omega=params->omega;
  double x0=params->x0;
  
  return -A*cos(omega*(x-x0))/omega;
}

int main()
{
  int workspace_size=1000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  double A=0.432,omega=0.7232,x0=0.1241;
  params_t params(A,omega,x0);
  
  //! function structure
  gsl_function f;
  f.function=kernel;
  f.params=&params;
  
  //integrate
  int method=1;
  double result;
  double abserr;
  double start=0.2345246,end=0.93245234,epsabs=0,epsrel=1e-6;
  
  //refernces: https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
  gsl_integration_qag(&f,start,end,epsabs,epsrel,workspace_size,method,workspace,&result,&abserr);
  
  //verify
  double analytic_result=kernel_primitive(end,&params)-kernel_primitive(start,&params);
  printf("Result: %.16lg , analytic: %.16lg, relative difference: %.16lg\n",result,analytic_result,result/analytic_result-1);
  
  gsl_integration_workspace_free(workspace);
  
  return 0;
}
