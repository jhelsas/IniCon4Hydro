#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "../splitandfit.h"

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
 
double gauss_e2s(double x[],size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }  
  
  double xx = x[0], yy = x[1],r2,r2hs;

  conv_wrap *cw = (conv_wrap *)(lpar->p);
  double *p = (double*)(cw->funpar);
  double    e0 = p[0];
  double   x0  = p[1];
  double   y0  = p[2];
  double  sigx = p[3];
  double  sigy = p[4];
  double  e0hs = p[5];
  double  x0hs = p[6];
  double  y0hs = p[7];
  double sigxhs = p[8];
  double sigyhs = p[9];
  double    tau = p[10];
  double epsilon,s;
  
  r2 = pow((xx-x0)/sigx,2.)+pow((yy-y0)/sigy,2.);
  r2hs= ((xx-x0hs)/sigxhs)*((xx-x0hs)/sigxhs) + ((yy-y0hs)/sigyhs)*((yy-y0hs)/sigyhs);

  if(r2hs > 20.)
    epsilon= e0*gsl_sf_exp(-0.5*r2);
    
  else
    epsilon= e0*gsl_sf_exp(-0.5*r2)+e0hs*gsl_sf_exp(-0.5*r2hs);
    
  s=cw->f2s(epsilon,cw->f2spar);
  return s*tau;
}

int null_velocity(double *x,size_t dim,void *par,double *u){
  int l;
  for(l=0;l<dim;l+=1)
    u[l]=0.;
  
  return 0;
}

int main(){
  int D=2,Ntri=6,split_type=0;
  int l,err,Npoints,N;
  double cutoff=0.02,xi[D],xf[D],p[11],xv[Ntri*(D+1)*D];
  double xl[D],xu[D],dx[D];
  double *xp,*x,*u,*S,s,dist,h=0.1;
  conv_wrap wrp;
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[0] =1.0*22;   /* e0 */   p[5] =0.0*22;  /* e0hs */
  p[1] =0.0;   /* x0 */   p[6] =0.7;  /* x0hs */
  p[2] =0.0;   /* y0 */   p[7] =1.25; /* y0hs */
  p[3] =1.0;  /* sigx */  p[8] =0.3; /* sigxhs */
  p[4] =1.5;  /* sigy */  p[9] =0.3; /* sigyhs */
  p[10] = 1.0; /* tau */
  
  F.f= gauss_e2s; F.dim=D;
  wrp.f2s=e2s_qg; wrp.f2spar=NULL; 
  wrp.funpar=(void*)p; par.p=(void*)&wrp;
  F.params=(void*)&par;
  
  if(split_type==0){
    cout << "init\n";
    for(l=0;l<D;l+=1){xi[l]=-4.0;xf[l]=4.0;}
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
  err=print_moving_sph(D,"results/gauss_e2s.dat",dom,null_velocity,&par); if(err!=0) return err;
  
  cout << "reading\n";
  err=sph_read("results/gauss_e2s.dat",&D,&N,&x,&u,&S);if(err!=0) return err;

  for(l=0;l<D;l+=1){xl[l]=-4.0;dx[l]=0.15;xu[l]=4.0+1.01*dx[l];}
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  
  cout << "ploting\n";
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,gauss_e2s,"results/gauss_e2s_plot.dat",&wrp);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S; 
  
  return 0;
}
