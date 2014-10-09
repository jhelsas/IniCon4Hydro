#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <vector>
#include "../splitandfit.h"
#include "../trial-functions.h"

using namespace std;

/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */

double gubser_entropy(double *x,size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2\n");
    abort();
  }
  double *p = (double*)(lpar->p);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  return (s0*(2.*q)*(2.*q)*(1.+q*q*(tau*tau+r2)))/(lambda*sqrt(lambda));
  
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

int main(){
  int D=2,Ntri=6,split_type=0;
  int l,err,Npoints,N;
  double cutoff=0.0002,xi[D],xf[D],p[3],xv[Ntri*(D+1)*D];
  double xl[D],xu[D],dx[D];
  double *xp,*x,*u,*S,s,dist,h=0.1;
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[0]=1.0; /*  s0 */ 
  p[1]=1.0; /*  q  */
  p[2]=1.0; /* tau */ 
  F.f= gubser_entropy; F.dim=D;par.p=(void*)p;F.params=(void*)&par;
  
  if(split_type==0){
    cout << "init\n";
    for(l=0;l<D;l+=1){xi[l]=-6.;xf[l]=6.;}
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
  err=print_moving_sph(D,"results/gubser.dat",dom,gubser_velocity,&par); if(err!=0) return err;
  
  cout << "reading\n";
  err=sph_read("results/gubser.dat",&D,&N,&x,&u,&S);if(err!=0) return err;

  for(l=0;l<D;l+=1){xl[l]=-6.0;dx[l]=0.15;xu[l]=6.0+1.01*dx[l];}
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  
  cout << "ploting\n";
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,gubser_entropy,"results/gubser_plot.dat",p);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S; 
  
  return 0;
}
