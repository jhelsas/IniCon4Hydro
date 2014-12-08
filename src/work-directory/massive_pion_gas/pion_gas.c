#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_integration.h>

double pKernel(double x,void *p){
  double mbT = *(double*)p;
  return x*sqrt(x*x-mbT*mbT)*log(1.-exp(-x));
}

double pressure(double mbT,void *params){
  size_t limit=1000,err;
  const double g=3;
  double *p,res,abserr;
  double eabs,erel,I,S,q=1.0,Q;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (10000);
  gsl_function F;
    
  eabs=0.0001;
  erel=1e-8;
  limit =10000;
  
  F.function = &pKernel;
  F.params = &mbT;
  err=gsl_integration_qagiu(&F,mbT,eabs,erel,limit,w,&I,&abserr);

  gsl_integration_workspace_free (w);
  
  return (-g*I)/(2*M_PI*M_PI);
}

double eKernel(double x,void *p){
  double mbT = *(double*)p;
  return (x*x*sqrt(x*x-mbT*mbT))/(exp(x)-1.);
}

double energy_density(double mbT,void *params){
  size_t limit=1000,err;
  const double g=3;
  double *p,res,abserr;
  double eabs,erel,I,S,q=1.0,Q;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (10000);
  gsl_function F;
    
  eabs=0.0001;
  erel=1e-8;
  limit =10000;
  
  F.function = &eKernel;
  F.params = &mbT;
  err=gsl_integration_qagiu(&F,mbT,eabs,erel,limit,w,&I,&abserr);

  gsl_integration_workspace_free (w);
  
  return (g*I)/(2*M_PI*M_PI);
}

int main(){
  int err;
  double m_pi=0.135,T,Ti=0.01,Tf=0.1,dT=0.01;
  double edens,press,entrop,soundsp,abserr;
  gsl_function Edens,Press;
  
  Edens.function = &energy_density;
  Edens.params = NULL;
  Press.function = &pressure;
  Press.params = NULL;
  
  for(T=Ti;T<=Tf;T+=dT){
    edens = energy_density(m_pi/T,NULL);
    press = pressure(m_pi/T,NULL);
    gsl_deriv_central(&Press,m_pi/T,1e-8,&entrop,&abserr);
    gsl_deriv_central(&Edens,m_pi/T,1e-8,&soundsp,&abserr);
    soundsp = entrop/soundsp;
    entrop += 3.0*press;
    printf("%lf %lf %lf %lf %lf\n",T,press,edens,entrop,soundsp);
  }
  
  return 0;
}
