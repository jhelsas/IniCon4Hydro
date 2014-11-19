#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct eosp{
  double T,s,p,e,h,hsh,cs2;
} eosp;

int eospS_landau(eosp *par){  
  const double C_pi = 0.134828384; /* 3*(hbarc)*((45÷(3x128×π^2))^(1/3)) GeV fm */
  par->p=C_pi*pow(par->s,4.0/3.0);
  par->e=3.0*(par->p);
  par->h=(par->e)+(par->p);
  par->hsh=(4.0/3.0)*(par->p);
  par->T=((par->h)/(par->s)); 
  return 0;
}

double e2s_pion(double epsilon,void *p){
  const double C_pi = 0.134828384; /* 3*(hbarc)*((45÷(3x128×π^2))^(1/3)) GeV fm */
  return pow(epsilon/(3.0*C_pi),3./4.);
}

int eospS_qg(eosp *par){  
  const double C_qg = 0.058356312; /* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */ 
  par->p=C_qg*pow(par->s,4.0/3.0);
  par->e=3.0*(par->p);
  par->h=(par->e)+(par->p);
  par->hsh=(4.0/3.0)*(par->p);
  par->T=((par->h)/(par->s)); 
  return 0;
}

double e2s_qg(double epsilon,void *p){
  const double C_qg = 0.058356312;/* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */ 
  return pow(epsilon/(3.0*C_qg),3./4.);
}

int eospS_qgphr(eosp *par){
  const double s1=2.1, s2=9.4125, C_hrg=0.1149;
  const double B=0.32, e1=0.28, e2=1.45;
  const double beta0=0.2, p1=0.056,Tc=(e1+p1)/s1;
  const double C_qgp = (e2-B)/pow(s2,4./3.);
  
  if(par->s < 0.)
    return 1;
  
  if( par->s < s1 ){
    par->e   = C_hrg*pow(par->s,1.+beta0);
    par->p   = C_hrg*beta0*pow(par->s,1.+beta0);
    par->T     = C_hrg*(1.+beta0)*pow(par->s,beta0);
    par->h   = C_hrg*(1.+beta0)*pow(par->s,1.+beta0);
    par->hsh = C_hrg*(1.+beta0)*beta0*pow(par->s,1.+beta0);
  }
  else if( s1<= par->s && par->s < s2 ){
    par->e   = Tc*(par->s)-p1;
    par->p   = p1;
    par->T     = (e1+p1)/s1;
    par->h   = (Tc/s1)*(par->s);
    par->hsh = 0.;
  }
  else {
    par->e   = C_qgp*pow(par->s,4./3.)+B;
    par->p   = (1./3.)*C_qgp*pow(par->s,4./3.)-B;
    par->T   = (4./3.)*C_qgp*pow(par->s,1./3.);
    par->h   = (4./3.)*C_qgp*pow(par->s,4./3.);
    par->hsh = (4./9.)*C_qgp*pow(par->s,4./3.);
  }
    
  return 0;
}

double e2s_qgphr(double epsilon, void *p){
  const double s1=2.1, s2=9.4125, C_hrg=0.1149;
  const double B=0.32, e1=0.28, e2=1.45;
  const double beta0=0.2, p1=0.056,Tc=(e1+p1)/s1;
  const double C_qgp = (e2-B)/pow(s2,4./3.);
  
  if(epsilon <= e1 )
    return pow(epsilon/C_hrg,1./(1.+beta0));
  else if(e1 < epsilon && epsilon <= e2 )
    return ((epsilon+p1)/(e1+p1))*s1;
  else
    return pow((epsilon-B)/(e2-B),3./4.)*s2;
}

int load_eos_table(char *filename,double *eos_t){
  const int Np=33001,Ne=5;
  int i;
  FILE *eosfile;
  double T,e,p,cs2,logs;
  eosfile = fopen(filename,"r");
  
  eos_t = (double*)malloc(Ne*Np*sizeof(double));
  
  if(eosfile==NULL)
    return 1;
    
  for(i=0;i<Np;i+=1){
    fscanf(eosfile,"%lf %lf %lf %lf %lf",&T,&cs2,&e,&p,&logs);
    eos_t[i*Ne+0] = logs; eos_t[i*Ne+1] = T;
    eos_t[i*Ne+2] = p;    eos_t[i*Ne+3] = e;
    eos_t[i*Ne+4] = cs2;
  }
  
  fclose(eosfile);
  return 0;
}

int EoS_table(eosp *par, double *eos_t){
  int i,Ne=5,Np=33001;
  double ds,dT,cs2,logs,sp,hp,pp,ep,T,dwds;
  double hbarc=0.197;
  
  ds = 0.0005;
  
  logs = log(par->s);
  if(logs < eos_t[0] || logs > eos_t[Ne*Np-1])
    return 1;
  i=(int)( (logs - eos_t[0])/ds );
  par->T = (eos_t[(i+1)*Ne+1]-eos_t[i*Ne+1])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+1];
  par->p = (eos_t[(i+1)*Ne+2]-eos_t[i*Ne+2])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+2];
  par->e = (eos_t[(i+1)*Ne+3]-eos_t[i*Ne+3])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+3];
    cs2  = (eos_t[(i+1)*Ne+4]-eos_t[i*Ne+4])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+4];
  par->h = par->e + par->p;
  
  dwds = (cs2+1.)*T;
  par->hsh=(4.0/3.0)*(par->p); /* verificar isso */
  
  return 0;
}

int main(){
  const double ei=8.14451245649177e-7,ef=322.538746367931;
  double logs,lsi,lsf,dls,*eos_t;
  eosp point;
  FILE *dadosout;
  
  lsi=log(e2s_qgphr(ei,NULL));
  lsf=log(e2s_qgphr(ef,NULL));
  
  dls=0.0005; // ds = 0.0005*hbarc;
  dadosout=fopen("qgphr_table.eos","w");  
  
  for(logs = lsi ; logs < lsf ; logs += dls){
    point.s = exp(logs);
    eospS_qgphr(&point);
    fprintf(dadosout,"%.10lf %.10lf %.10lf %.10lf %.10lf\n",point.T,point.cs2,point.e,point.p,logs);
  }
  fclose(dadosout);
  
  if(eos_t!=NULL)
    free(eos_t);
  
  return 0;
}
