#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

double f(double *x,size_t n,void *params){
  if(fabs(x[0]) <= 1. )
    return 0.75*(1.-(x[0]*x[0]));
  else
    return 0;
}

double fd(double *x,size_t n,void *params){
  unsigned int i;
  double r=0.;
  for(i=0;i<=n;i+=1)
    r+=x[i]*x[i];
  r=sqrt(r);
  if(r<=1. && n==1)
    return 0.75*(1.-r*r);
  else if(r<=1. && n==2)
    return (3.*M_1_PI)*(1.-r*r);
  else if(r<=1. && n==3)
    return (1.875*M_1_PI)*(1.-r*r);
  else
    return 0;
}

double gd(double *x,size_t n,void *params){
  unsigned int i;
  double re=1.;
  for(i=0;i<n;i+=1)
    if(fabs(x[i])>1.)
      return 0;
  for(i=0;i<n;i+=1)
    re *= 0.75*(1-x[i]*x[i]);
  return re;
}

double inicon(double x[],size_t dim,void * par) {

  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }  

  double xx = x[0];
  double yy = x[1];

  double *p = (double*)par;
  double    A = p[0];
  double   x0 = p[1];
  double   y0 = p[2];
  double sigx = p[3];
  double sigy = p[4];

  return A*exp(-0.5*(pow((xx-x0)/sigx,2.)+pow((yy-y0)/sigy,2.)));

}
