#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

typedef struct EoSstable{
  int Ne,Np;
  double *eos_t;
}  EoSstable;

typedef struct eosp{
  double T,s,p,e,h,hsh,cs2;
} eosp;

typedef struct SPHeq_particle
{
  int D,id,*ind,fo; 
  double ni,rho,q,S,s,s_p,sigma,Nb,Nc;
  double rho_p,e_p,p_p,h_p,hsh_p,T,nb,nc,cs2;
  double *x,*u,*v,*grads;
  double Ta,rhopa,spa,*xa,*ua,*gradsa;
  double Tb,rhopb,spb,*xb,*ub,*gradsb;
  double Tc,rhopc,spc,*xc,*uc,*gradsc;
} SPHeq_particle;

int eospS_landau(eosp *par,void* params){  
  const double C_pi = 0.134828384; 
  par->p=C_pi*pow(par->s,4.0/3.0);
  par->e=3.0*(par->p);
  par->h=(par->e)+(par->p);
  par->hsh=(4.0/3.0)*(par->p);
  par->T=((par->h)/(par->s)); 
  par->cs2 = 1./3.;
  return 0;
}

int EoS_pi(SPHeq_particle *par)
{	  
  const double Cpi = 0.134828384; 
  const double etacharg_qg;
  par->p_p=Cpi*pow(par->s_p,4.0/3.0);
  par->e_p=3.0*(par->p_p);
  par->h_p=(par->e_p)+(par->p_p);
  par->hsh_p=(4.0/3.0)*(par->p_p);
  par->T=((par->h_p)/(par->s_p)); 
  par->cs2 = 1./3.;
  par->nc = (par->q)*etacharg_qg*(par->s_p); 
  par->Nc = (par->nc)/(par->rho);
  
  return 0;
}

double e2s_pion(double epsilon,void *p){
  const double C_pi = 0.134828384; /* 3*(hbarc)*((45÷(3x128×π^2))^(1/3)) GeV fm */
  return pow(epsilon/(3.0*C_pi),3./4.);
}

int eospS_qg(eosp *par,void* params){  
  const double C_qg = 0.0179276286445;
  //const double C_qg = 0.01948439; /* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */ 
  par->p=C_qg*pow(par->s,4.0/3.0);
  par->e=3.0*(par->p);
  par->h=(par->e)+(par->p);
  par->hsh=(4.0/3.0)*(par->p);
  par->T=((par->h)/(par->s)); 
  par->cs2 = 1./3.;
  return 0;
}

int EoS_qg(SPHeq_particle *par)
{	  
  //const double Cqg = 0.01948439;  
  const double Cqg = 0.0179276286445;
  const double etacharg_qg;
  par->p_p=Cqg*pow(par->s_p,4.0/3.0);
  par->e_p=3.0*(par->p_p);
  par->h_p=(par->e_p)+(par->p_p);
  par->hsh_p=(4.0/3.0)*(par->p_p);
  par->T=((par->h_p)/(par->s_p)); 
  par->cs2 = 1./3.;
  par->nc = (par->q)*etacharg_qg*(par->s_p); 
  par->Nc = (par->nc)/(par->rho);
  
  return 0;
}

double e2s_qg(double epsilon,void *p){
  //const double C_qg = 0.01948439; /* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */  
  const double C_qg = 0.0179276286445;
  return pow(epsilon/(3.0*C_qg),3./4.);
}

int eospS_qgphr(eosp *par,void* params){
  const double s1=2.1, s2=9.4125, C_hrg=0.1149; // (s1,s2) alternativo : (2.4,10.7571)
  const double B=0.32, e1=0.28, e2=1.45;
  const double beta0=0.2, p1=0.056,Tc=(e1+p1)/s1; // Tc atual: 0.160 
  const double C_qgp = (e2-B)/pow(s2,4./3.); // livro : 0.0538961130616 | atual: 0.0568597733934
  
  if(par->s < 0.)
    return 1;
  
  if( par->s < s1 ){
    par->e   = C_hrg*pow(par->s,1.+beta0);
    par->p   = C_hrg*beta0*pow(par->s,1.+beta0);
    par->T     = C_hrg*(1.+beta0)*pow(par->s,beta0);
    par->h   = C_hrg*(1.+beta0)*pow(par->s,1.+beta0);
    par->hsh = C_hrg*(1.+beta0)*beta0*pow(par->s,1.+beta0);
    par->cs2 = (par->hsh)/(par->h);
  }
  else if( s1<= par->s && par->s < s2 ){
    par->e   = Tc*(par->s)-p1;
    par->p   = p1;
    par->T     = (e1+p1)/s1;
    par->h   = (Tc/s1)*(par->s);
    par->hsh = 0.;
    par->cs2 = (par->hsh)/(par->h);
  }
  else {
    par->e   = C_qgp*pow(par->s,4./3.)+B;
    par->p   = (1./3.)*C_qgp*pow(par->s,4./3.)-B;
    par->T   = (4./3.)*C_qgp*pow(par->s,1./3.);
    par->h   = (4./3.)*C_qgp*pow(par->s,4./3.);
    par->hsh = (4./9.)*C_qgp*pow(par->s,4./3.);
    par->cs2 = (par->hsh)/(par->h);
  }
    
  return 0;
}

int EoS_qgp(SPHeq_particle *par){
  const double s1=2.1, s2=9.4125, C_hrg=0.1149;
  const double B=0.32, e1=0.28, e2=1.45;
  const double beta0=0.2, p1=0.056,Tc=(e1+p1)/s1;
  const double C_qgp = (e2-B)/pow(s2,4./3.);
  
  if(par->s_p < 0.)
    return 1;
  
  if( par->s_p < s1 ){
    par->e_p   = C_hrg*pow(par->s_p,1.+beta0);
    par->p_p   = C_hrg*beta0*pow(par->s_p,1.+beta0);
    par->T     = C_hrg*(1.+beta0)*pow(par->s_p,beta0);
    par->h_p   = C_hrg*(1.+beta0)*pow(par->s_p,1.+beta0);
    par->hsh_p = C_hrg*(1.+beta0)*beta0*pow(par->s_p,1.+beta0);
    par->cs2 = (par->hsh_p)/(par->h_p);
  }
  else if( s1<= par->s_p && par->s_p < s2 ){
    par->e_p   = Tc*(par->s_p)-p1;
    par->p_p   = p1;
    par->T     = (e1+p1)/s1;
    par->h_p   = (Tc/s1)*(par->s_p);
    par->hsh_p = 0.;
    par->cs2 = (par->hsh_p)/(par->h_p);
  }
  else {
    par->e_p   = C_qgp*pow(par->s_p,4./3.)+B;
    par->p_p   = (1./3.)*C_qgp*pow(par->s_p,4./3.)-B;
    par->T     = (4./3.)*C_qgp*pow(par->s_p,1./3.);
    par->h_p   = (4./3.)*C_qgp*pow(par->s_p,4./3.);
    par->hsh_p = (4./9.)*C_qgp*pow(par->s_p,4./3.);
    par->cs2 = (par->hsh_p)/(par->h_p);
  }
  par->nc = 0.0; 
  par->Nc = 0.0;
    
  return 0;
}

double e2s_qgphr(double epsilon, void *p){
  const double s1=2.1, s2=9.4125, C_hrg=0.1149;
  const double B=0.32, e1=0.28, e2=1.45;
  const double beta0=0.2, p1=0.056,Tc=(e1+p1)/s1;
  const double C_qgp = (e2-B)/pow(s2,4./3.); // 0.0189532577978
  
  if(epsilon <= e1 )
    return pow(epsilon/C_hrg,1./(1.+beta0));
  else if(e1 < epsilon && epsilon <= e2 )
    return ((epsilon+p1)/(e1+p1))*s1;
  else
    return pow((epsilon-B)/(e2-B),3./4.)*s2;
}

int load_eos_table(char *filename,double **eos_t){
  const int Np=33001,Ne=5;
  int i;
  FILE *eosfile;
  double T,e,p,cs2,logs;
  eosfile = fopen(filename,"r");
  
  *eos_t = (double*)malloc(Ne*Np*sizeof(double));
  
  if(eosfile==NULL)
    return 1;
    
  for(i=0;i<Np;i+=1){
    fscanf(eosfile,"%lf %lf %lf %lf %lf",&T,&cs2,&e,&p,&logs);
    (*eos_t)[i*Ne+0] = logs; (*eos_t)[i*Ne+1] = T;
    (*eos_t)[i*Ne+2] = p;    (*eos_t)[i*Ne+3] = e;
    (*eos_t)[i*Ne+4] = cs2;
  }
  
  fclose(eosfile);
  return 0;
}

int EoS_table(eosp *par, double *eos_t){
  int i,Ne=5,Np=33001;
  double ds,dT,cs2,logs,sp,hp,pp,ep,T,dwds;
  double hbarc=0.1973269718;
  
  ds = 0.0005;
  
  logs = log(par->s);
  
  if(logs < eos_t[0] || logs > eos_t[Ne*(Np-1)]){
    printf("%lf %lf %lf\n", eos_t[0],logs,eos_t[Ne*(Np-1)]);
    return 1;}
    
  i=(int)( (logs - eos_t[0])/ds );
  par->T = (eos_t[(i+1)*Ne+1]-eos_t[i*Ne+1])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+1];
  par->p = (eos_t[(i+1)*Ne+2]-eos_t[i*Ne+2])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+2];
  par->e = (eos_t[(i+1)*Ne+3]-eos_t[i*Ne+3])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+3];
  par->cs2  = (eos_t[(i+1)*Ne+4]-eos_t[i*Ne+4])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+4];
  par->h = par->e + par->p;
  par->hsh=(par->cs2)*(par->h); 
  
  return 0;
}

int load_zoltan_table(char *filename,double *z_table){
  int i,Nl=18,Ne=5;
  double T,p,dp,emp,demp,cs2,dcs2,e,de,s,ds;
  const double hbarc=0.1973269718;
  FILE *dadosin;
  
  dadosin=fopen(filename,"r");
  if(dadosin==NULL)
    return 1;
  for(i=0;i<Nl;i+=1){
    fscanf(dadosin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&T,&p,&dp,&emp,&demp,&cs2,&dcs2,&e,&de,&s,&ds);
    T=T/1000.;
    p = (p*T*T*T*T)/(hbarc*hbarc*hbarc);
    e = (e*T*T*T*T)/(hbarc*hbarc*hbarc);
    s = (s*T*T*T)/(hbarc*hbarc*hbarc);
    
    z_table[i*Ne+0] = log(s);   z_table[i*Ne+1] = T;
    z_table[i*Ne+2] = p;        z_table[i*Ne+3] = e;
    z_table[i*Ne+4] = cs2;
  }
  fclose(dadosin);
  return 0;
}

double lag2w(double x,double xa,double xb,double xc){
  return ((x-xa)*(x-xb))/((xc-xa)*(xc-xb));
}
  
int zoltan_eos_old(eosp *par, double *eos_t){
  int good=1,i,iup,idown,Ne=5,Np=18,i1,i2,i3;
  double logs,dwds;
  double ls1,ls2,ls3,f1,f2,f3;
  double hbarc=0.1973269718;
  
  logs = log(par->s);
  if(logs < eos_t[0] || logs > eos_t[Ne*(Np-1)] )
    return 1;
  
  idown=0;
  iup=Np-1;
  while(good!=0){
    i = (iup+idown)/2;
    if(logs < eos_t[i*Ne])
      iup  = i;
    else if(logs > eos_t[(i+1)*Ne])
      idown = i;
    else{
      if(i == 0){
        ls1 = eos_t[Ne*0]; i1 = 0;
        ls2 = eos_t[Ne*1]; i2 = 1;
        ls3 = eos_t[Ne*2]; i3 = 2;
      } else if( i == (Np-2) ){
        ls1 = eos_t[Ne*(Np-3)]; i1 = Np-3;
        ls2 = eos_t[Ne*(Np-2)]; i2 = Np-2;
        ls3 = eos_t[Ne*(Np-1)]; i3 = Np-1;
      } else {
        ls1 = eos_t[Ne*(i-1)]; i1 = i-1; 
        ls2 = eos_t[Ne*i];     i2 = i;
        ls3 = eos_t[Ne*(i+1)]; i3 = i+1;
      }  
      good = 0;
    }
  }
  f1 = lag2w(logs,ls2,ls3,ls1);
  f2 = lag2w(logs,ls3,ls1,ls2);
  f3 = lag2w(logs,ls1,ls2,ls3);
  
  par->T = f1*eos_t[i1*Ne+1] + f2*eos_t[i2*Ne+1] + f3*eos_t[i3*Ne+1];
  par->p = f1*eos_t[i1*Ne+2] + f2*eos_t[i2*Ne+2] + f3*eos_t[i3*Ne+2];
  par->e = f1*eos_t[i1*Ne+3] + f2*eos_t[i2*Ne+3] + f3*eos_t[i3*Ne+3];
  par->cs2=f1*eos_t[i1*Ne+4] + f2*eos_t[i2*Ne+4] + f3*eos_t[i3*Ne+4];
  par->h = par->e + par->p; 
  par->hsh=(par->cs2)*(par->h); 
    
  return 0;
}

gsl_interp_accel *acc; // global variables
gsl_spline *splT, *splp ,*sple ,*splcs2 ;

int eos_zoltan_old2(eosp *par, void *params){
  static int call_count = 0;
  const int Ne=5,Np=18;
  const double hbarc=0.1973269718;
  const double ls[18] = {-1.769049, -1.077466, -0.395522, 0.142758, 
	                      0.549412, 0.824491, 1.108151, 1.425858, 
	                      1.732869, 2.388880, 2.923345, 3.268573, 
	                      3.909100, 4.582473, 5.569299, 6.140774, 
	                      7.031874, 7.722272 };
  const double T[18]  = {0.100000, 0.115000, 0.129000, 0.139000, 
	                     0.147000, 0.152000, 0.158000, 0.166000, 
	                     0.175000, 0.200000, 0.228000, 0.250000, 
	                     0.299000, 0.366000, 0.500000, 0.600000, 
	                     0.800000, 1.000000 };
  const double e[18]  = {0.014186, 0.032551, 0.073524, 0.137981, 
	                     0.221213, 0.302902, 0.419333, 0.602841, 
	                     0.858120, 1.844991, 3.499477, 5.333056, 
	                     11.848111, 27.884914, 100.540059, 212.359345, 
	                     686.086922, 1707.554136 };
  const double p[18]  =  {0.002863, 0.006601, 0.013335, 0.022349, 
	                      0.033425, 0.043768, 0.059210, 0.087956, 
	                      0.131831, 0.335264, 0.742100, 1.235398, 
	                      3.058248, 7.893719, 30.585002, 66.288501, 
	                      219.633110, 550.530030 };
  const double cs2[18]=  {0.190000, 0.180000, 0.140000, 0.130000, 
	                      0.120000, 0.120000, 0.140000, 0.160000, 
	                      0.180000, 0.220000, 0.260000, 0.270000, 
	                      0.290000, 0.320000, 0.320000, 0.320000, 
	                      0.320000, 0.320000 };
  double logs;
  
  if(call_count==0){
    acc=gsl_interp_accel_alloc ();
    //splT = gsl_spline_alloc (gsl_interp_steffen, 18);
    //splp = gsl_spline_alloc (gsl_interp_steffen, 18);
    //sple = gsl_spline_alloc (gsl_interp_steffen, 18);
    //splcs2 = gsl_spline_alloc (gsl_interp_steffen, 18);
    splT = gsl_spline_alloc (gsl_interp_cspline, 18);
    splp = gsl_spline_alloc (gsl_interp_cspline, 18);
    sple = gsl_spline_alloc (gsl_interp_cspline, 18);
    splcs2 = gsl_spline_alloc (gsl_interp_cspline, 18);
    gsl_spline_init (splT, ls, T, Np);
    gsl_spline_init (splp, ls, p, Np);
    gsl_spline_init (sple, ls, e, Np);
    gsl_spline_init (splcs2, ls, cs2, Np);  
    call_count+=1;
  }
  
  if(par->s<0)
    return (-1);
    
  logs = log(par->s);
  if(logs < ls[0] || logs > ls[Np-1])
    return 1;
  
  par->T = gsl_spline_eval (splT, logs, acc);
  par->p = gsl_spline_eval (splp, logs, acc);
  par->e = gsl_spline_eval (sple, logs, acc);
  par->cs2 = gsl_spline_eval (splcs2, logs, acc);
  par->h = par->e + par->p; 
  par->hsh=(par->cs2)*(par->h); 
  
  return 0;
}

int eos_zoltan(eosp *par, void *params){
  static int call_count = 0;
  const int Ne=5,Np=18;
  const double hbarc=0.1973269718;
  const double ls[18] = {-1.769049, -1.077466, -0.395522, 0.142758, 
	                      0.549412, 0.824491, 1.108151, 1.425858, 
	                      1.732869, 2.388880, 2.923345, 3.268573, 
	                      3.909100, 4.582473, 5.569299, 6.140774, 
	                      7.031874, 7.722272 };
	                      
  const double T[18]  = {0.100, 0.115, 0.129, 0.139,0.147, 0.152, 0.158, 
	                     0.166, 0.175, 0.200, 0.228,0.250, 0.299, 0.366, 
	                     0.500, 0.600,0.800, 1.000 };
  const double e[18]  = {1.09,1.43,2.04,2.84,3.64,4.36,5.17,6.1,7.03,
	                     8.86,9.95,10.49,11.39,11.94,12.36,12.59,
	                     12.87,13.12 };
  const double p[18]  =  {0.22,0.29,0.37,0.46,0.55,0.63,0.73,0.89,
	                      1.08,1.61,2.11,2.43,2.94,3.38,3.76,3.93,
	                      4.12,4.23 };
  const double cs2[18]=  {0.19,0.18,0.14,0.13,0.12,0.12,0.14,0.16,
	                      0.18,0.22,0.26,0.27,0.29,0.32,0.32,0.32,
	                      0.32,0.32 };
  double logs;
  
  if(call_count==0){
    acc=gsl_interp_accel_alloc ();
    //splT = gsl_spline_alloc (gsl_interp_steffen, 18);
    //splp = gsl_spline_alloc (gsl_interp_steffen, 18);
    //sple = gsl_spline_alloc (gsl_interp_steffen, 18);
    //splcs2 = gsl_spline_alloc (gsl_interp_steffen, 18);
    splT = gsl_spline_alloc (gsl_interp_cspline, 18);
    splp = gsl_spline_alloc (gsl_interp_cspline, 18);
    sple = gsl_spline_alloc (gsl_interp_cspline, 18);
    splcs2 = gsl_spline_alloc (gsl_interp_cspline, 18);
    gsl_spline_init (splT, ls, T, Np);
    gsl_spline_init (splp, ls, p, Np);
    gsl_spline_init (sple, ls, e, Np);
    gsl_spline_init (splcs2, ls, cs2, Np);  
    call_count+=1;
  }
  
  if(par->s<0)
    return (-1);
    
  logs = log(par->s);
  if(logs < ls[0] || logs > ls[Np-1])
    return 1;
  
  par->T = gsl_spline_eval (splT, logs, acc);
  par->p = gsl_spline_eval (splp, logs, acc);
  par->e = gsl_spline_eval (sple, logs, acc);
  par->cs2 = gsl_spline_eval (splcs2, logs, acc);
  par->p *= (par->T)*(par->T/hbarc)*(par->T/hbarc)*(par->T/hbarc);
  par->e *= (par->T)*(par->T/hbarc)*(par->T/hbarc)*(par->T/hbarc);
  par->h = par->e + par->p; 
  par->hsh=(par->cs2)*(par->h); 
  
  return 0;
}

int EoS_zoltan(SPHeq_particle *par){
  static int call_count = 0;
  const int Ne=5,Np=18;
  const double hbarc=0.1973269718;
  const double ls[18] = {-1.769049, -1.077466, -0.395522, 0.142758, 
	                      0.549412, 0.824491, 1.108151, 1.425858, 
	                      1.732869, 2.388880, 2.923345, 3.268573, 
	                      3.909100, 4.582473, 5.569299, 6.140774, 
	                      7.031874, 7.722272 };
	                      
  const double T[18]  = {0.100, 0.115, 0.129, 0.139,0.147, 0.152, 0.158, 
	                     0.166, 0.175, 0.200, 0.228,0.250, 0.299, 0.366, 
	                     0.500, 0.600,0.800, 1.000 };
  const double e[18]  = {1.09,1.43,2.04,2.84,3.64,4.36,5.17,6.1,7.03,
	                     8.86,9.95,10.49,11.39,11.94,12.36,12.59,
	                     12.87,13.12 };
  const double p[18]  =  {0.22,0.29,0.37,0.46,0.55,0.63,0.73,0.89,
	                      1.08,1.61,2.11,2.43,2.94,3.38,3.76,3.93,
	                      4.12,4.23 };
  const double cs2[18]=  {0.19,0.18,0.14,0.13,0.12,0.12,0.14,0.16,
	                      0.18,0.22,0.26,0.27,0.29,0.32,0.32,0.32,
	                      0.32,0.32 };
  double logs;
  
  if(call_count==0){
    acc=gsl_interp_accel_alloc ();
    //splT = gsl_spline_alloc (gsl_interp_steffen, 18);
    //splp = gsl_spline_alloc (gsl_interp_steffen, 18);
    //sple = gsl_spline_alloc (gsl_interp_steffen, 18);
    //splcs2 = gsl_spline_alloc (gsl_interp_steffen, 18);
    splT = gsl_spline_alloc (gsl_interp_cspline, 18);
    splp = gsl_spline_alloc (gsl_interp_cspline, 18);
    sple = gsl_spline_alloc (gsl_interp_cspline, 18);
    splcs2 = gsl_spline_alloc (gsl_interp_cspline, 18);
    gsl_spline_init (splT, ls, T, Np);
    gsl_spline_init (splp, ls, p, Np);
    gsl_spline_init (sple, ls, e, Np);
    gsl_spline_init (splcs2, ls, cs2, Np);  
    call_count+=1;
  }
  
  if(par->s_p<0)
    return (-1);
    
  logs = log(par->s_p);
  if(logs < ls[0] || logs > ls[Np-1])
    return 1;
    
  
  par->T = gsl_spline_eval (splT, logs, acc);
  par->p_p = gsl_spline_eval (splp, logs, acc);
  par->e_p = gsl_spline_eval (sple, logs, acc);
  par->cs2 = gsl_spline_eval (splcs2, logs, acc);
  par->p_p *= (par->T)*(par->T/hbarc)*(par->T/hbarc)*(par->T/hbarc);
  par->e_p *= (par->T)*(par->T/hbarc)*(par->T/hbarc)*(par->T/hbarc);
  par->h_p = par->e_p + par->p_p; 
  par->hsh_p=(par->cs2)*(par->h_p); 
  
  par->nc = 0.0; 
  par->Nc = 0.0;
    
  return 0;
}

double e2s_zoltanRF(double s,void *params){
  int err;
  double e0=*(double*)params;
  eosp par;
  par.s=s;
  err = eos_zoltan(&par,NULL);
  if(err!=0)
    return (-1);
  
  return (par.e - e0);
}

const gsl_root_fsolver_type *rsT;
gsl_root_fsolver *solver=NULL;
    
double e2s_zoltan_old(double e, void *p){
  static int call_count = 0;
  int status, iter, max_iter=100;
  double lsi = exp(-1.769049)-400.,lsf=exp(6.5)+400., ei=0.014186, ef=1707.554136;
  double x_lo,x_hi;
  double delta,a=-7.73192e-07,b=0.000927288,c=0.194341,d=-0.164113;
  gsl_function F;
  double s;
    
  if(call_count==0){
    rsT = gsl_root_fsolver_brent;
    solver = gsl_root_fsolver_alloc (rsT);
    call_count +=1;
  }
  /*
  if(e<5.)
    delta=20.;
  else if(e < 10.)
    delta=35.;
  else if(e < 20.)
    delta=50.;
  else if(e<50.)
    delta= 100.;
  else if(e<100.)
    delta = 200.;
  else 
    delta = 400;
  
  lsi = a*(e*e*e)+b*(e*e)+c*e+e - delta;
  if(lsi < exp(-1.769049))
    lsi=0.;
  lsf = a*(e*e*e)+b*(e*e)+c*e+e + delta;*/
    
  F.function = &e2s_zoltanRF;
  F.params = &e;
  gsl_root_fsolver_set (solver, &F, lsi, lsf);
  do{
    iter++;
    status = gsl_root_fsolver_iterate (solver);
    s = gsl_root_fsolver_root (solver);
    x_lo = gsl_root_fsolver_x_lower (solver);
    x_hi = gsl_root_fsolver_x_upper (solver);
    status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
    if (status == GSL_SUCCESS)
      break;
  } while (status == GSL_CONTINUE && iter < max_iter);
  if(iter>=max_iter)
    printf("max iterations reached\n");  
  
  return s;
}

gsl_interp_accel *global_e2s_zoltan_acc; // global variables
gsl_spline *global_e2s_zoltan_spls;

double e2s_zoltan(double e,void *ep){
  const int Np=83;
  static int static_call_count=0;
  double ez_table[83] = {0.01,0.26,0.51,0.76,1.01,1.26,1.51,1.76,2.01,
                         2.26,2.51,2.76,3.01,3.26,3.51,3.76,4.01,4.26,
                         4.51,4.76,5.00,7.00,9.00,11.00,13.00,15.00,
                         17.0,19.0,21.0,23.0,25.0,27.0,29.0,31.0,33.0,
                         35.0,37.0,39.0,41.0,43.0,45.00,52.5,60.0,
                         67.5,75.0,82.5,90.0,97.5,105.0,112.50,120.0,
                         127.5,135.0,142.5,150.0,157.5,165.0,172.5,
                         180.0,187.5,195.0,202.5,210.0,217.5,225.0,
                         232.5,240.0,247.5,255.0,262.5,270.0,277.5,
                         285.0,292.5,300.0,307.5,315.0,322.5,330.0,
                         337.5,345.0,352.5,360.0};
  double sz_table[83] = {0.170524,1.994627,3.594898,5.091375,6.513195,
                         7.878641,9.196814,10.474862,11.718528,
                         12.932077,14.118870,15.281767,16.423192,
                         17.545205,18.649566,19.737736,20.810858,
                         21.869890,22.915678,23.948969,24.929789,
                         32.741990,40.052551,46.997441,53.661896,
                         60.097748,66.338096,72.408171,78.327858,
                         84.113186,89.777319,95.331205,100.784158,
                         106.144529,111.419782,116.616567,121.740296,
                         126.796059,131.788246,136.720760,141.597103,
                         159.433296,176.653015,193.352457,209.603689,
                         225.462452,240.973086,256.171676,271.087773,
                         285.747860,300.175215,314.389804,328.408896,
                         342.247557,355.919018,369.434987,382.805883,
                         396.041033,409.148822,422.136823,435.011905,
                         447.780318,460.447772,473.019361,485.498318,
                         497.886828,510.186980,522.400789,534.530511,
                         546.577365,558.543514,570.430711,582.240652,
                         593.974976,605.635268,617.223062,628.739843,
                         640.187046,651.566062,662.878239,674.124883,
                         685.307257,696.426589};
  double s;
  
  if(e < ez_table[0] || e > ez_table[Np-1])
    return (-1.0);  
  
  if(static_call_count==0){
    global_e2s_zoltan_acc=gsl_interp_accel_alloc ();
    global_e2s_zoltan_spls = gsl_spline_alloc (gsl_interp_cspline, 83);
    gsl_spline_init (global_e2s_zoltan_spls, ez_table, sz_table, 83);
    static_call_count+=1;
  }
  
  s = gsl_spline_eval (global_e2s_zoltan_spls, e, global_e2s_zoltan_acc);
  
  return s;
}
 
int print_zoltan(){
  int err;
  double lsi,s2,lsf,logs,dls,z_table[5*18];
  double e,de,ei=0.01,ef=360;
  const double hbarc=0.1973269718;
  SPHeq_particle par; 
  eosp point;
  FILE *dadosout;
  /*
  lsi=-1.769049;//log(e2s_qgphr(ei,NULL));
  lsf= 6.5;     //log(e2s_qgphr(ef,NULL));
    
  err= load_zoltan_table("zoltan.dat",z_table);
  if(err!=0){printf("zoltan table not loaded\n");return err;}
  dls=0.0005; 
  
  
  dadosout=fopen("results/zoltan_table.eos","w");  
  for(logs = lsi ; logs < lsf ; logs += dls){
    point.s = exp(logs);
    par.s_p = point.s;
    err=eos_zoltan(&point,z_table);
    err=EoS_zoltan(&par);
    if(err!=0)
      continue;
    fprintf(dadosout,"%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf ",
            point.T,point.cs2,point.e,point.p,logs,
            (pow(hbarc,3)*point.e)/pow(point.T,4),
            (pow(hbarc,3)*point.p)/pow(point.T,4));
    fprintf(dadosout,"%lf %lf %lf %lf %lf %lf %lf\n",par.T,par.cs2,par.e_p,par.p_p,
            (pow(hbarc,3)*par.e_p)/pow(par.T,4),
            (pow(hbarc,3)*par.p_p)/pow(par.T,4));
  }
  fclose(dadosout); */
  
  de=0.01;
  dadosout = fopen("results/zoltan_e2s_comparison.dat","w");
  for(e=ei;e<=ef;e+=de){
    logs = e2s_zoltan_old(e,NULL); 
    s2 = e2s_zoltan(e,NULL);
    point.s = logs;
    err=eos_zoltan(&point,NULL);
    if(err!=0)
      continue;
    fprintf(dadosout,"%lf %lf %lf %lf\n",e,logs,s2,logs-s2);
    //fprintf(dadosout,"%lf %lf \n",e,s2);
  }
  fclose(dadosout);
  
  /*
    
  de=0.25;
  ei=0.01; ef=5.0;
  dadosout = fopen("results/zoltan_e2s_array.dat","w");
  for(e=ei;e<=ef;e+=de){
    logs = e2s_zoltan(e,NULL);
    point.s = logs;
    err=eos_zoltan(&point,NULL); 
    if(err!=0)
      continue;
    fprintf(dadosout,"%lf %lf \n",e,logs);
  }
  de=2.0;
  ei=ef; ef = 45.;
  for(e=ei;e<=ef;e+=de){
    logs = e2s_zoltan(e,NULL);
    point.s = logs;
    err=eos_zoltan(&point,NULL); 
    if(err!=0)
      continue;
    fprintf(dadosout,"%lf %lf \n",e,logs);
  }
  de=7.5;
  ei=ef+de; ef = 360.;
  for(e=ei;e<=ef;e+=de){
    logs = e2s_zoltan(e,NULL);
    point.s = logs;
    err=eos_zoltan(&point,NULL); 
    if(err!=0)
      continue;
    fprintf(dadosout,"%lf %lf \n",e,logs);
  }
  fclose(dadosout);
  
  gsl_spline_free(splT); 
  gsl_spline_free(sple); 
  gsl_spline_free(splp); 
  gsl_spline_free(splcs2);
  gsl_interp_accel_free(acc); 
  */
  
  if(solver!=NULL)
    gsl_root_fsolver_free (solver);
    
  return 0;
} 
 
int print_pasi(){
  int err;
  double lsi,lsf,logs,dls,*eos_t;
  double e,de,ei=0.01,ef=360;
  const double hbarc=0.1973269718;
  eosp point;
  SPHeq_particle par; 
  FILE *dadosout;
  
  lsi=-1.769049;//log(e2s_qgphr(ei,NULL));
  lsf= 6.5;     //log(e2s_qgphr(ef,NULL));
  dls=0.0005; // ds = 0.0005*hbarc;
  
  err=load_eos_table("EoS_pasi.eos",&eos_t);
  if(err!=0){printf("not loading pasi EoS\n");return 1;}
  
  dadosout=fopen("results/pasi_table.eos","w");  
  for(logs = lsi ; logs < lsf ; logs += dls){
    point.s = exp(logs);
    par.s_p = point.s;
    err=EoS_table(&point,eos_t); 
    if(err!=0){printf("skiped %lf\n",logs);continue;}
    fprintf(dadosout,"%.10lf %.10lf %.10lf %.10lf %.10lf\n",point.T,point.cs2,point.e,point.p,logs);
  }
  fclose(dadosout); 
  
  return 0;
}
 
int print_qg(){
  int err;
  double lsi,lsf,logs,dls,z_table[5*18];
  double e,de,ei=0.01,ef=360;
  const double hbarc=0.1973269718;
  eosp point;
  SPHeq_particle par; 
  FILE *dadosout;
  
  lsi=-1.769049;//log(e2s_qgphr(ei,NULL));
  lsf= 6.5;     //log(e2s_qgphr(ef,NULL));
  dls=0.0005; // ds = 0.0005*hbarc;
    
  dadosout=fopen("results/qg_table.eos","w");  
  for(logs = lsi ; logs < lsf ; logs += dls){
    point.s = exp(logs);
    err=eospS_qg(&point,z_table); 
    err=EoS_qg(&par);
    //if(err!=0)
    //  continue;
    fprintf(dadosout,"%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf ",
            point.T,point.cs2,point.e,point.p,logs,
            (pow(hbarc,3)*point.e)/pow(point.T,4),
            (pow(hbarc,3)*point.p)/pow(point.T,4));
    fprintf(dadosout,"%lf %lf %lf %lf %lf %lf %lf\n",par.T,par.cs2,par.e_p,par.p_p,
            (pow(hbarc,3)*par.e_p)/pow(par.T,4),
            (pow(hbarc,3)*par.p_p)/pow(par.T,4));
  }
  fclose(dadosout); 
  
  de=0.01;
  dadosout = fopen("results/qg_e2s_comparison.dat","w");
  for(e=ei;e<=ef;e+=de){
    logs = e2s_qg(e,NULL);
    point.s = logs;
    err=eospS_qg(&point,NULL); 
    //if(err!=0)
    //  continue;
    fprintf(dadosout,"%lf %lf %lf %lf %lf\n",logs,e,point.e,e-point.e,point.T);
  }
  fclose(dadosout);
  
  return 0;
} 
 
int print_qgphr(){
  int err;
  double lsi,lsf,logs,dls,z_table[5*18];
  double e,de,ei=0.01,ef=360;
  const double hbarc=0.1973269718;
  eosp point;
  SPHeq_particle par; 
  FILE *dadosout;
  
  lsi=-1.769049;//log(e2s_qgphr(ei,NULL));
  lsf= 6.5;     //log(e2s_qgphr(ef,NULL));
  dls=0.0005; // ds = 0.0005*hbarc;
    
  dadosout=fopen("results/qgphr_table.eos","w");  
  for(logs = lsi ; logs < lsf ; logs += dls){
    
    point.s = exp(logs);
    par.s_p = point.s;
    err=eospS_qgphr(&point,NULL); 
    err=EoS_qgp(&par);
    if(err!=0)
      continue;
    fprintf(dadosout,"%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf ",
            point.T,point.cs2,point.e,point.p,logs,
            (pow(hbarc,3)*point.e)/pow(point.T,4),
            (pow(hbarc,3)*point.p)/pow(point.T,4));
    fprintf(dadosout,"%lf %lf %lf %lf %lf %lf %lf\n",par.T,par.cs2,par.e_p,par.p_p,
            (pow(hbarc,3)*par.e_p)/pow(par.T,4),
            (pow(hbarc,3)*par.p_p)/pow(par.T,4));
  }
  fclose(dadosout); 
  
  de=0.01;
  dadosout = fopen("results/qgphr_e2s_comparison.dat","w");
  for(e=ei;e<=ef;e+=de){
    logs = e2s_qgphr(e,NULL);
    point.s = logs;
    err=eospS_qgphr(&point,NULL); 
    if(err!=0)
      continue;
    fprintf(dadosout,"%lf %lf %lf %lf %lf\n",logs,e,point.e,e-point.e,point.T);
  }
  fclose(dadosout);
  
  return 0;
} 

int main(){
  int i,Np=18,Ne=5,err;
  const double ei=0.01,ef=360.;
  double e,de=0.005,logs,lsi,lsf,dls,*eos_t=NULL,z_table[5*18];
  eosp point;
  FILE *dadosout;
  
  err=print_zoltan();
  //err=print_pasi();
  //err=print_qg();
  //err=print_qgphr();
  
  if(eos_t!=NULL)
    free(eos_t);
  
  return 0;
}
