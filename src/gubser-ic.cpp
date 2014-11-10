#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

double gubser_energy2entropy(double *x,size_t dim, void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2\n");
    abort();
  }
  conv_wrap *cw = (double *)(lpar->p);
  double *p = (double*)(cw->funpar);
  double e0=p[0], q=p[1],tau=p[2];
  double r2,lambda,gamma,s,epsilon;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  epsilon = (e0*pow(2.*q,8./3.))/(tau*lambda*cbrt(tau*lambda));
  s=cw->f2s(epsilon,cw->f2spar);
  gamma = (1.+q*q*(tau*tau+r2))/sqrt(lambda);
  return s*gamma*tau;
}

int gubser_velocity(double *x,size_t dim,void *par,double *u){
  if(dim!=2)
    return 1;
    
  int err;
  wparams *lpar=(wparams*)par;
  double *p = (double*)(lpar->p);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda,c;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  c=(2.0*q*q*tau)/sqrt(lambda);
  u[0]=c*x[0];
  u[1]=c*x[1];
  
  return 0;
}
