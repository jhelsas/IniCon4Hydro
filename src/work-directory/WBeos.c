#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_deriv.h>

/*
 * BW stands for Budapest-Wuppertal Collaboration
 * 
 * parametrization can be found at 
 *    arXiv:1007.2580v2 [hep-lat]
 * 
 * and updated parameters at
 *    arXiv: 1309.5258 [hep-lat]
 * 
 * functions below calculate:
 * 
 *  trace defect : (e-3p)/T^4
 *    pressure   :   p/T^4
 * energy density:   e/T^4
 * entropy dens. :   s/T^3
 *  sound speed  :   c_s^2
 * 
 *   in admensional values. 
 */

int BW_traceDefect(double T,void *params){
  const double h0=0.1396,h1=-0.18,h2=0.035;
  const double f0=1.05,f1=6.39,f2=-4.72;
  const double g1=-0.92,g2=0.57;
  double I=0.,t=(T/0.2); /* t= T/200 MeV*/
  
  I = gsl_sf_lnsinh(f1*t+f2)-gsl_sf_lncosh(f1*t+f2);
  I = f0*(gsl_sf_exp(I)+1.);
  I = I/(1.+g1*t+g2*t*t);
  I *= gsl_sf_exp(-h1/t-h2/(t*t));
  
  return I;
}

double BW_dIdT(double T,void *params){
  const double h0=0.1396,h1=-0.18,h2=0.035;
  const double f0=1.05,f1=6.39,f2=-4.72;
  const double g1=-0.92,g2=0.57;
  double dIdt=0.,F0,t=(T/0.2);
  
  dIdt = h0+(f0*tanh(f1*t+f2))/(1.+g1*t+g2*t*t);
  dIdt *= h1/t+h2/(t*t);
  dIdt += (f0*f1*sech(f1*t+f2)*sech(f1*t+f2));
  dIdt -= (f0*(tanh(f1*t+f2)+1.)*(g1+g2*t))/(1+g1*t+g2*t*t);
  dIdt *= exp(-h1/t-h2/(t*t));
  
  return dIdt/(0.2*0.2*0.2*0.2*0.2);  
}

int BW_pressure(double T,void *params){
  static gsl_integration_workspace* w=NULL ;
  static int callcount = 0;
  gsl_function F;
  int err;
  double pressure=0.,abserr;
  
  if(callcount <=0){
    w=gsl_integration_workspace(1000);
    F.function=&BW_traceDefect;
    F.params=NULL;
  }
  if(T<=0 && w !=NULL){
    gsl_integration_workspace_free(w);
    return -1.;
  }
  
  gsl_integration_qag(&F,0.,T.,1e-8,1000,w,&pressure,&err);
  
  return pressure;
}

int BW_energyDen(double T,void *params){
  static gsl_integration_workspace* w=NULL ;
  static int callcount = 0;
  gsl_function F;
  int err;
  double pressure=0.,abserr,I;
  
  if(callcount <=0){
    w=gsl_integration_workspace(1000);
    F.function=&BW_traceDefect;
    F.params=NULL;
  }
  if(T<=0 && w !=NULL){
    gsl_integration_workspace_free(w);
    return -1.;
  }
  
  I=BW_traceDefect(T,NULL);
  gsl_integration_qag(&F,0.,T.,1e-8,1000,w,&pressure,&err);
  
  return 3.0*pressure+I;
}

int BW_entropyDen(double T,void *params){
  static gsl_integration_workspace* w=NULL ;
  static int callcount = 0;
  gsl_function F;
  int err;
  double pressure=0.,energyDens,abserr,I;
  
  if(callcount <=0){
    w=gsl_integration_workspace(1000);
    F.function=&BW_traceDefect;
    F.params=NULL;
  }
  if(T<=0 && w !=NULL){
    gsl_integration_workspace_free(w);
    return -1.;
  }
  
  I=BW_traceDefect(T,NULL);
  gsl_integration_qag(&F,0.,T.,1e-8,1000,w,&pressure,&err);
  energyDens=3.0*pressure+I;  
  
  return energyDens+pressure;
}

int main(){
  int i;
  double T,Tmin,Tmax,dT;
  double eT4,pT4,sT3, cs2;
  
  
  return 0;
}
