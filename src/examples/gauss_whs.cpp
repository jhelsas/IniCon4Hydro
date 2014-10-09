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
 
double gauss_whs(double x[],size_t dim,void *par){
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

  double *p = (double*)(lpar->p);
  double    A = p[0];
  double   x0 = p[1];
  double   y0 = p[2];
  double sigx = p[3];
  double sigy = p[4];
  double    Ahs = p[5];
  double   x0hs = p[6];
  double   y0hs = p[7];
  double sigxhs = p[8];
  double sigyhs = p[9];
  
  r2 = pow((xx-x0)/sigx,2.)+pow((yy-y0)/sigy,2.);
  r2hs= ((xx-x0hs)/sigxhs)*((xx-x0hs)/sigxhs) + ((yy-y0hs)/sigyhs)*((yy-y0hs)/sigyhs);

  if(r2hs > 20.)
    return A*gsl_sf_exp(-0.5*r2);
    
  else
    return A*gsl_sf_exp(-0.5*r2)+Ahs*gsl_sf_exp(-0.5*r2hs);
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
  double cutoff=0.0002,xi[D],xf[D],p[10],xv[Ntri*(D+1)*D];
  double xl[D],xu[D],dx[D];
  double *xp,*x,*u,*S,s,dist,h=0.1;
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  p[0] =1.0;p[5] =0.7;
  p[1] =0.0;p[6] =0.7;
  p[2] =0.0;p[7] =1.25;
  p[3] =1.0;p[8] =0.3;
  p[4] =1.5;p[9] =0.3;
  
  F.f= gauss_whs; 
  F.dim=D;
  par.p=(void*)p;
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
  err=print_moving_sph(D,"results/gauss_whs.dat",dom,null_velocity,&par); if(err!=0) return err;
  
  cout << "reading\n";
  err=sph_read("results/gauss_whs.dat",&D,&N,&x,&u,&S);if(err!=0) return err;

  for(l=0;l<D;l+=1){xl[l]=-4.0;dx[l]=0.15;xu[l]=4.0+1.01*dx[l];}
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  
  cout << "ploting\n";
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,gauss_whs,"results/gauss_whs_plot.dat",p);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S; 
  
  return 0;
}
