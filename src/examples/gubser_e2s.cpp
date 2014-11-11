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
typedef struct conv_wrap{
  double (*f2s)(double,void*p);
  void *funpar, *f2spar;
} conv_wrap;

double e2s_qg(double epsilon,void *p){
  const double C_qg = 0.058356312;/* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */ 
  return pow(epsilon/C_qg,3./4.);
}

double gubser_e2s(double *x,size_t dim, void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2\n");
    abort();
  }
  conv_wrap *cw = (conv_wrap *)(lpar->p);
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
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2\n");
    abort();
  }
  conv_wrap *cw = (conv_wrap *)(lpar->p);
  double *p = (double*)(cw->funpar);
  double e0=p[0], q=p[1],tau=p[2];
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
  double cutoff,xi[D],xf[D],p[3],xv[Ntri*(D+1)*D];
  double xl[D],xu[D],dx[D];
  double *xp,*x,*u,*S,s,dist,h=0.1;
  conv_wrap wrp;
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[0]=50.689088426;/*  e0 */ 
  p[1]=1.0; /*  q  */
  p[2]=1.0; /* tau */ 
  F.f= gubser_e2s; F.dim=D;
  wrp.f2s=e2s_qg; wrp.f2spar=NULL; 
  wrp.funpar=(void*)p; par.p=(void*)&wrp;
  F.params=(void*)&par;
  cutoff=0.03202;
  
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
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,gubser_e2s,"results/gubser_plot.dat",&wrp);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S; 
  
  return 0;
}
