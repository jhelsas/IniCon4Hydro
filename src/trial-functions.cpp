#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "trial-functions.h"
#include "splitandfit.h"

using namespace std;

/*
 * Exemplos de função a ser utilizada com o código inicon4hydro
 * 
 * Toda função para ser usada como input do código, precisa ter o
 * formato igual às funções exemplo abaixo:
 * 
 * double func(double *x,size_t dim,void *par){
 *   [suas variaveis]
 *   int err;
 *   wparams *lpar=(wparams*)par;
 * 
 *   err=check_inside(x,n,lpar->mdel);
 *   if(err!=0)
 *     return 0.;
 * 
 *   [resto do teu código normal]
 * }
 *  
 * É absolutamente fundamental que o código tenha o formato acima 
 * prescrito, por questão de compatibilidade com a GSL e com o 
 * algorítimo de particionamento de domínios do código. É necessário
 * também utilizar wparams como parâmetro de entrada, colocando em seu
 * membro p o resto dos parâmetros que se desejar usar, no mesmo jeito
 * que se faria numa função habitual da GSL.
 * 
 * Eu cheguei a escrever um wrapper para que fosse apenas 
 * necessário manter o formato de tipos da GSL, e internamente o código
 * de particionamento de domínios se encarregasse de checar se o ponto
 * em questão está ou não dentro do domínio, mas isso teve uma 
 * penalidade grande de performance, portanto, por enquanto, requisito
 * que seja seguido o modelo acima.
 */

double fd(double *x,size_t n,void *par){
  unsigned int i;
  int err;
  double r=0.;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,n,lpar->mdel);
  if(err!=0)
    return 0.;
    
  for(i=0;i<=n;i+=1)
    r+=x[i]*x[i];
  r=sqrt(r);
  if(r<=1. && n==1)
    return 0.75*(1.-r*r);
  else if(r<=1. && n==2)
    return (3.*M_1_PI)*(1.-r*r);
  else if(r<=1. && n==3)
    return (1.875*M_1_PI)*(1.-r*r);
  else
    return 0;
  
}

double winicon(double x[],size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }  
  
  double xx = x[0];
  double yy = x[1];

  double *p = (double*)(lpar->p);
  double    A = p[0];
  double   x0 = p[1];
  double   y0 = p[2];
  double sigx = p[3];
  double sigy = p[4];

  return A*exp(-0.5*(pow((xx-x0)/sigx,2.)+pow((yy-y0)/sigy,2.)));
}

double woodsaxon(double x[],size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }  
  
  int l;
  double r=0.,wds=0.;
  double *p = (double*)(lpar->p);
  double    A = p[0];
  double   x0 = p[1];
  double   y0 = p[2];
  double sigx = p[3];
  double sigy = p[4];
  double    a = p[5];
  double    R = p[6];
  
  r=sqrt( (sigx*(x[0]-x0))*(sigx*(x[0]-x0)) + (sigy*(x[1]-y0))*(sigy*(x[1]-y0)) ) ;
  return A/(1.0+exp(a*(r-R)));
}

double woodsaxon_whotspot(double x[],size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }  
  
  int l;
  double r=0.,r2hs=0.;
  double *p = (double*)(lpar->p);
    
  double    Aws = p[0];double    Ahs = p[7];
  double   x0ws = p[1];double   x0hs = p[8];
  double   y0ws = p[2];double   y0hs = p[9];
  double sigxws = p[3];double sigxhs = p[10];
  double sigyws = p[4];double sigyhs = p[11];
  double    aws = p[5];
  double    Rws = p[6];
   
  r=sqrt( ((x[0]-x0ws)/sigxws)*((x[0]-x0ws)/sigxws) + ((x[1]-y0ws)/sigyws)*((x[1]-y0ws)/sigyws) ) ;
  r2hs= ((x[0]-x0hs)/sigxhs)*((x[0]-x0hs)/sigxhs) + ((x[1]-y0hs)/sigyhs)*((x[1]-y0hs)/sigyhs);
  
  return Aws/(1.0+gsl_sf_exp((r-Rws)/aws)) + Ahs*gsl_sf_exp(-0.5*r2hs);
}

double gubser_entropy(double x[],size_t dim,void *par){
  int err;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,dim,lpar->mdel);
  if(err!=0)
    return 0.;
  
  if(dim != 2) {
    fprintf(stderr, "error: dim != 2");
    abort();
  }
  double *p = (double*)(lpar->p);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda;
  r2=x[0]*x[0]+x[1]*x[1];
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  return (s0*(2.*q)*(2.*q))/(tau*lambda);
  
}

int gubser_velocity(double x[],size_t dim,void *par,double *u){
  if(dim!=2)
    return 1;
  
  double *p = (double*)(lpar->p);
  double s0=p[0], q=p[1],tau=p[2];
  double r2,lambda;
  lambda = 1.+2.*q*q*(tau*tau+r2) + q*q*q*q*(tau*tau-r2)*(tau*tau-r2);
  
  if(r2<=0){
    u[0]=0.;u[1]=0.;
    return 0;
  }
  
  u[0]=x[0]/(sqrt(lambda*r2));
  u[1]=x[1]/(sqrt(lambda*r2));
  
  return 0;
}
