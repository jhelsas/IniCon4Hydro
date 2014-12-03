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
  const double C_qg = 0.0179276286445;/* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */ 
  return pow(epsilon/(3.0*C_qg),3./4.);
}

double rot_conf_soliton_e2s(double *x,size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 3) {
    fprintf(stderr, "error: dim != 3\n");
    abort();
  }
  conv_wrap *cw = (conv_wrap *)(lpar->p);
  double *p = (double*)(cw->funpar);
  double e0 = p[0], w = p[1], t0=p[2], L = p[3];
  double r,rT2,lambda,gamma,s,epsilon;
  rT2 = x[0]*x[0] + x[1]*x[1];
  r = sqrt(rT2 + x[2]*x[2]);
  lambda = (L*L+(t0+r)*(t0+r))*(L*L+(t0-r)*(t0-r)) - 4.*w*w*L*L*rT2;
  gamma = (L*L+r*r+t0*t0)/sqrt(lambda);
  epsilon = e0/(lambda*lambda);
  
  s=cw->f2s(epsilon,cw->f2spar);
  return s*gamma;
}

int rot_conf_soliton_vel(double *x,size_t dim, void *par,double *u){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return err;
  
  if(dim != 3) {
    fprintf(stderr, "error: dim != 3\n");
    abort();
  }
  conv_wrap *cw = (conv_wrap *)(lpar->p);
  double *p = (double*)(cw->funpar);
  double e0 = p[0], w = p[1], t0=p[2], L = p[3];
  double r,rT2,lambda,gamma,s,epsilon;
  rT2 = x[0]*x[0] + x[1]*x[1];
  r = sqrt(rT2 + x[2]*x[2]);
  lambda = (L*L+(t0+r)*(t0+r))*(L*L+(t0-r)*(t0-r)) - 4.*w*w*L*L*rT2;
  
  u[0] = ( 2.*t0*x[0] + 2.*w*L*x[1] )/sqrt(lambda);
  u[1] = ( 2.*t0*x[1] - 2.*w*L*x[0] )/sqrt(lambda);
  u[2] = ( 2.*t0*x[2] )/sqrt(lambda);
  
  return 0.;
}

int main(){
  int D=3,Ntri=6,split_type=0;
  int l,err,Npoints,N;
  double cutoff,xi[D],xf[D],p[4],xv[Ntri*(D+1)*D];
  double xl[D],xu[D],dx[D];
  double *xp,*x,*u,*S,s,dist,h=0.1;
  conv_wrap wrp;
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[1] = 0.20; /* w obs: 0 <= w <= 1 */
  p[2] = 0.00; /* t0 */
  p[3] = 2.50; /* L */
  p[0] = 22.0*pow(p[3],8); /* e0 obs: >=0 */
  
  F.f= rot_conf_soliton_e2s; F.dim=D;
  wrp.f2s=e2s_qg; wrp.f2spar=NULL; 
  wrp.funpar=(void*)p; par.p=(void*)&wrp;
  F.params=(void*)&par;
  cutoff=0.125;
  
  if(split_type==0){
    cout << "init\n";
    for(l=0;l<D;l+=1){xi[l]=-8.;xf[l]=8.;}
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
  
  //cout << "clean\n";
  //err=clean_domain(dom);if(err!=0) return err;
  
  cout << "print\n";  
  err=print_moving_sph(D,"results/rotconfsol/rotconfsol.dat",dom,rot_conf_soliton_vel,&par); if(err!=0) return err;
  
  cout << "reading\n";
  err=sph_read("results/rotconfsol/rotconfsol.dat",&D,&N,&x,&u,&S);if(err!=0) return err;

  for(l=0;l<D;l+=1){xl[l]=-8.0;dx[l]=0.2;xu[l]=8.0+1.01*dx[l];}
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  
  cout << "ploting\n";
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,rot_conf_soliton_e2s,"results/rotconfsol/rotconfsol_plot.dat",&wrp);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S; 
  
  return 0;
}

