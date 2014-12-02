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

  //TH3D *hTable = (TH3D*)(lpar->p);
  TH3D **p = (TH3D**)(lpar->p);
  TH3D *hEdens = *&p[0];
  TH3D *hGamma = *&p[1];
  double zmin = hEdens->GetZaxis()->GetBinCenter(1);
  double zmax = hEdens->GetZaxis()->GetBinCenter(hEdens->GetNbinsZ());
  double edens = 0.;
  double gamma = 1.;
  if(zz>zmin && zz<zmax) {
    edens = (double)(hEdens->Interpolate(xx,yy,zz));
    gamma = (double)(hGamma->Interpolate(xx,yy,zz));
  }
  //printf("x: %2.8lf   y: %2.8lf   z: %2.8lf     e: %2.8lf     g: %2.8lf\n",xx,yy,zz,edens,gamma);

  //value = cw->f2s(value,cw->f2spar);
  return gamma*edens;

}

int phsd_velocity(double *x,size_t dim,void *par,double *u){
  if(dim!=3)
    return 1;
    
  double xx = x[0];
  double yy = x[1];
  double zz = x[2];

  int err;
  wparams *lpar=(wparams*)par;
  TH3D **p = (TH3D**)(lpar->p);
  TH3D *v = 0;

  for(int i=0;i<3;++i) {
    v = *&p[i];
    u[i] = (double)v->Interpolate(xx,yy,zz);
  }
  
  return 0;
}

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
    wparams vpar;
    vector <domain> dom;
    gsl_monte_function F;
    ifstream sphfile;
    ofstream plotfile;
    FILE *sphofile;

    // get phsd snapshot "tables" (TH3D ROOT histograms)
    TFile *fICo = new TFile("phsd-ico_NUM1_t0.12.root","READ");
    TH3D *p[2];
    p[0] = (TH3D*)fICo->Get("hEdensXYZ");
    p[1] = (TH3D*)fICo->Get("hGammaXYZ");
    par.p=(void*)p; 

    TH3D *v[3];
    v[0] = (TH3D*)fICo->Get("hBetaX");
    v[1] = (TH3D*)fICo->Get("hBetaY");
    v[2] = (TH3D*)fICo->Get("hBetaZ");
    vpar.p=(void*)v;
  
    F.f= phsd_edens; 
    F.dim=D;
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
    //err=print_moving_sph(D,"results/phsd_ico.dat",dom,null_velocity,&par); if(err!=0) return err;
    err=print_moving_sph(D,"results/phsd_ico.dat",dom,phsd_velocity,&vpar); if(err!=0) return err;

    cout << "reading\n";
    err=sph_read("results/phsd_ico.dat",&D,&N,&x,&u,&S);if(err!=0) return err;

    for(l=0;l<2;l+=1) { xl[l]=-10.0; dx[l]=0.2; xu[l]=10.0+1.01*dx[l]; }

    err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         

    cout << "ploting\n";
    err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,phsd_edens,"results/phsd_ico_plot.dat",&p); if(err!=0){ cout << err << endl; return err;}

    //hEdens = 0;
    fICo->Close();

    delete x; delete u; delete S; 

    return 0;
}

