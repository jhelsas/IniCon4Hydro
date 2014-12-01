#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "../splitandfit.h"

// ROOT libraries
#include "TFile.h"
#include "TH3D.h"
#include "TApplication.h"

// this function returns the energy density at a given point (x,y,z)
// by interpolating the points of a table (TH3D) passed as parameter
double phsd_edens(double *x,size_t dim, void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 3) {
    fprintf(stderr, "error: dim != 3\n");
    abort();
  }

  double xx = x[0];
  double yy = x[1];
  double zz = x[2];

  TH3D *hTable = (TH3D*)(lpar->p);
  double zmin = hTable->GetZaxis()->GetBinCenter(1);
  double zmax = hTable->GetZaxis()->GetBinCenter(hTable->GetNbinsZ());
  double value = 0.;
  if(zz>zmin && zz<zmax) value = (double)(hTable->Interpolate(xx,yy,zz));
  //printf("x: %2.8lf   y: %2.8lf   z: %2.8lf     e: %2.8lf\n",xx,yy,zz,value);

  //value = cw->f2s(value,cw->f2spar);
  return value;

}

//int phsd_velocity(double *x,size_t dim,void *par,double *u){
//  if(dim!=2)
//    return 1;
//    
//  int err;
//  wparams *lpar=(wparams*)par;
//  double *p = (double*)(lpar->p);
//  double s0=p[0], q=p[1],tau=p[2];
//  double r2,lambda,c;
//  r2=x[0]*x[0]+x[1]*x[1];
//  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
//  c=(2.0*q*q*tau)/sqrt(lambda);
//  u[0]=c*x[0];
//  u[1]=c*x[1];
//  
//  return 0;
//}

int null_velocity(double *x,size_t dim,void *par,double *u){
    int l;
    for(l=0;l<dim;l+=1)
        u[l]=0.;

    return 0;
}


int main(){

    // this is needed to compile using ROOT libs
    TApplication *a = new TApplication("a",0,0);

    // define function properties
    int D=3,Ntri=6,split_type=0;
    int l,err,Npoints,N;
    double cutoff=2,xi[D],xf[D],xv[Ntri*(D+1)*D];
    double xl[D],xu[D],dx[D];
    double *xp,*x,*u,*S,s,dist,h=0.1;
    wparams par;
    vector <domain> dom;
    gsl_monte_function F;
    ifstream sphfile;
    ofstream plotfile;
    FILE *sphofile;

    // get phsd snapshot "table" (TH3D ROOT histogram)
    TFile *fICo = new TFile("phsd-ico_NUM30_t013.root","READ");
    TH3D *hEdens = (TH3D*)fICo->Get("hEdensXYZ");
  
    F.f= phsd_edens; 
    F.dim=D;
    par.p=(void*)hEdens; 
    F.params=(void*)&par;

    if(split_type==0){
        cout << "init\n";
        for(l=0;l<D;l+=1){xi[l]=-10.0;xf[l]=10.0;}
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
    err=print_moving_sph(D,"results/phsd_ico.dat",dom,null_velocity,&par); if(err!=0) return err;

    cout << "reading\n";
    err=sph_read("results/phsd_ico.dat",&D,&N,&x,&u,&S);if(err!=0) return err;

    for(l=0;l<2;l+=1) { xl[l]=-10.0; dx[l]=0.2; xu[l]=10.0+1.01*dx[l]; }

    err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         

    cout << "ploting\n";
    err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,phsd_edens,"results/phsd_ico_plot.dat",hEdens); if(err!=0){ cout << err << endl; return err;}

    hEdens = 0;
    fICo->Close();

    delete x; delete u; delete S; 

    return 0;
}

