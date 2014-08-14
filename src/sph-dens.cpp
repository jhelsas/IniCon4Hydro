#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

#define _USE_MATH_DEFINES

#define DIM 2

//#define A_d 1.  
 #define A_d (15.)/(7.*M_PI) 
// #define A_d (3.0)/(2.0*MYPI) 

double inicon(double x[],size_t dim,void * par) {

  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }  

  double xx = x[0];
  double yy = x[1];

  double *p = (double*)par;
  double    A = p[0];
  double   x0 = p[1];
  double   y0 = p[2];
  double sigx = p[3];
  double sigy = p[4];

  return A*exp(-0.5*(pow((xx-x0)/sigx,2.)+pow((yy-y0)/sigy,2.)));

}

double w_bspline(double r,double h)
{
	double R;
	if(h<=0.) exit(10);
	R=fabs(r)/h;
	if(R>=2.)
		return 0;
	else if((1.<=R)&&(R<2.))
		return (A_d)*(1./6.)*(2.-R)*(2.-R)*(2.-R)/pow(h,DIM);
	else
		return ((A_d)*((2./3.)-(R*R) + (R*R*R/2.0)))/pow(h,DIM);
}

int main(){
  const int N=2608,D=2,Npoints=161*161;
  int i,j,k,l;
  double p[5];
  double xp[Npoints*D],x[N*D],u[N*D],hp[N],rho[N],S[N],s,X[D],dist,h,a,b,c;
  ifstream sphfile,gridfile;
  ofstream plotfile,splashfile;
  
  p[0]=1.; 
  p[1]=0.;
  p[2]=0.;
  p[3]=.5;
  p[4]=.5;
  
  h=0.15;
  
  sphfile.open("SPH-particles.dat");
  for(i=0;i<N;i+=1){
    for(l=0;l<D;l+=1)
      sphfile >> x[i*D+l];
    for(l=0;l<D;l+=1)
      sphfile >> u[i*D+l];
    sphfile >> S[i];
  }
  gridfile.open("points-trans.gfcp");
  for(i=0;i<Npoints;i+=1){
    gridfile >> a >> b >> c;
    xp[D*i+0] = a/20.;
    xp[D*i+1] = b/20.;
  }
  gridfile.close();
    
  plotfile.open("ploting.dat");
  for(k=0;k<Npoints;k+=1){
    s=0.;
    for(i=0;i<N;i+=1){
      dist=0.;
      for(l=0;l<D;l+=1)
        dist+=(xp[k*D+l]-x[i*D+l])*(xp[k*D+l]-x[i*D+l]);
      dist=sqrt(dist);
      
      s+=S[i]*w_bspline(dist,h);
    }
    a=inicon(xp+k*D,2,p);
    for(l=0;l<D;l+=1)
      plotfile << xp[k*D+l] << " ";
    plotfile << s << " " << a << "\n";
  }
  plotfile.close();
  
  return 0;
}
