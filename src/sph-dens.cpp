#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include "splitandfit.h"
#include "trial-functions.h"

using namespace std;

#define _USE_MATH_DEFINES

int main(){
  int N=0,D=0;
  int i,j,k,l,err,Npoints;
  double p[12];
  double *xp,*x,*u,*S,s,dist,h,a,b,c;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[0]=1.; p[1]=0.; p[2]=0.;
  p[3]=1.; p[4]=1.; p[5]=1./25.; p[6]=1.;  
  p[7]=0.; p[8]=0.2; p[9]=0.2; p[10]=0.1; p[11]=0.1;
  
  h=0.05;    
    
  err=sph_read("SPH-particles.dat",&D,&N,&x,&u,&S);if(err!=0) return err;
  double xl[D],xu[D],dx[D];
  for(l=0;l<D;l+=1){xl[l]=-1.2;dx[l]=0.05;xu[l]=1.2+1.01*dx[l];}
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,woodsaxon_whotspot,"ploting.dat",p);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S;
  
  return 0;
}
