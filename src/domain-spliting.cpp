#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <vector>
#include "splitandfit.h"
#include "trial-functions.h"

using namespace std;

/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */

int main(int argc,char **argv){
  const int D=2,Ntri=6,split_type=0;
  int err,l; 
  double cutoff=0.0002,xi[D],xf[D],p[3],xv[Ntri*(D+1)*D];
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  
  p[0]=1.0; /*  s0 */ 
  p[1]=1.0; /*  q  */
  p[2]=1.0; /* tau */ 
  F.f= gubser_entropy; F.dim=D;par.p=(void*)p;F.params=(void*)&par;
  
  if(split_type==0){
    cout << "init\n";
    for(l=0;l<D;l+=1){xi[l]=-3.;xf[l]=3.;}
    err=init_cube(xi,xf,dom,D);if(err!=0) return err; 
  }
  else if(split_type==1){
    cout << "unit\n";
    err=unit2_hexagon(Ntri,D,xv);if(err!=0) return err;
    
    cout << "init\n";
    err=init_triangle(D,Ntri,xv,dom);if(err!=0) return err;
  } 
  
  cout << "split\n";
  err=domain_split(D,cutoff,dom,F); if(err!=0){ cout << "out: " << err << endl;return err;}
  
  cout << "clean\n";
  err=clean_domain(dom);if(err!=0) return err;
  
  cout << "print\n";  
  //err=print_sph(D,"SPH-particles.dat",dom); if(err!=0) return err;
  err=print_moving_sph(D,"SPH-particles.dat",dom,gubser_velocity,&par); if(err!=0) return err;
  
  return 0;
}
