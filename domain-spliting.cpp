#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "trial-functions.h"
#include "splitandfit.h"
#include <string.h>
#include <vector>

#define _USE_MATH_DEFINES

using namespace std;

/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */
 
int main(){
  const int D=2;
  int err,l,Nsplit=1<<D;
  double cutoff=0.001;
  double xi[D],xf[D],p[5];
  
  domain mdel;
  vector <domain> dom;
  gsl_monte_function F;
  
  p[0]=1.; p[1]=0.; p[2]=0.;
  p[3]=.5; p[4]=.5;
  F.f = &inicon; F.dim=D; F.params=(void*)p;
  
  for(l=0;l<D;l+=1){xi[l]=-1.;xf[l]=1.;}
  
  err=init_cube(xi,xf,dom,D);if(err!=0) return err;
  
  err=domain_split(D,Nsplit,cutoff,dom,F); if(err!=0) return err;
  
  err=clean_domain(dom);if(err!=0) return err;
  
  err=print_sph(D,"SPH-particles.dat",dom); if(err!=0) return err;
    
  return 0;
}
