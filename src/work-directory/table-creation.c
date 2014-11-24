#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

typedef struct EoSstable{
  int Ne,Np;
  double *eos_t;
}  EoSstable;

typedef struct eosp{
  double T,s,p,e,h,hsh,cs2;
} eosp;

int eospS_landau(eosp *par,void* params){  
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

int eospS_qg(eosp *par,void* params){  
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

int eospS_qgphr(eosp *par,void* params){
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
  double hbarc=0.1973269718;
  
  ds = 0.0005;
  
  logs = log(par->s);
  if(logs < eos_t[0] || logs > eos_t[Ne*(Np-1)])
    return 1;
  i=(int)( (logs - eos_t[0])/ds );
  par->T = (eos_t[(i+1)*Ne+1]-eos_t[i*Ne+1])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+1];
  par->p = (eos_t[(i+1)*Ne+2]-eos_t[i*Ne+2])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+2];
  par->e = (eos_t[(i+1)*Ne+3]-eos_t[i*Ne+3])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+3];
  par->cs2  = (eos_t[(i+1)*Ne+4]-eos_t[i*Ne+4])*( (logs-eos_t[i*Ne])/ds ) + eos_t[i*Ne+4];
  par->h = par->e + par->p;
  
  dwds = (par->cs2+1.)*(par->T);
  par->hsh=(4.0/3.0)*(par->p); /* verificar isso */
  
  return 0;
}



int load_zoltan_table(char *filename,double *z_table){
  int i,Nl=18,Ne=5;
  double T,p,dp,emp,demp,cs2,dcs2,e,de,s,ds;
  const double hbarc=0.1973269718;
  FILE *dadosin;
  
  dadosin=fopen(filename,"r");
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
  
int zoltan_eos(eosp *par, double *eos_t){
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
  
  dwds = (par->cs2+1.)*(par->T);
  par->hsh=(4.0/3.0)*(par->p); /* verificar isso */
    
  return 0;
}

int zoltan_eos2(eosp *par, double *eos_t){
  static int call_count=0;
  const int Ne=5,Np=18;
  int good=1,i,i1,i2,i3;
  double logs,dwds;
  double ls1,ls2,ls3,f1,f2,f3;
  double hbarc=0.1973269718;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *splT = gsl_spline_alloc (gsl_interp_cspline, Np);
  gsl_spline *splp = gsl_spline_alloc (gsl_interp_cspline, Np);
  gsl_spline *sple = gsl_spline_alloc (gsl_interp_cspline, Np);
  gsl_spline *splcs2 = gsl_spline_alloc (gsl_interp_cspline, Np);
  double T[Np],e[Np],p[Np],cs2[Np],ls[Np];
  
  for(i=0;i<Np;i+=1){
   ls[i] = eos_t[i*Ne]; T[i] = eos_t[i*Ne+1];
   p[i] = eos_t[i*Ne+2];e[i] = eos_t[i*Ne+3]; cs2[i] = eos_t[i*Ne+4];
  }
    
  gsl_spline_init (splT, ls, T, Np);
  gsl_spline_init (splp, ls, p, Np);
  gsl_spline_init (sple, ls, e, Np);
  gsl_spline_init (splcs2, ls, cs2, Np);
    
  logs = log(par->s);
  if(logs < eos_t[0] || logs > eos_t[Ne*(Np-1)] )
    return 1;
  
  par->T = gsl_spline_eval (splT, logs, acc);
  par->p = gsl_spline_eval (splp, logs, acc);
  par->e = gsl_spline_eval (sple, logs, acc);
  par->cs2 = gsl_spline_eval (splcs2, logs, acc);
  par->h = par->e + par->p; 
  
  dwds = (par->cs2+1.)*(par->T);
  par->hsh=(4.0/3.0)*(par->p); /* verificar isso */
  
  gsl_spline_free(splT); 
  gsl_spline_free(sple); 
  gsl_spline_free(splp); 
  gsl_spline_free(splcs2);
  gsl_interp_accel_free(acc); 
    
  call_count+=1;
    
  return 0;
}


#define EoS zoltan_eos
#define e2s e2s_zoltan

int main(){
  int i,Np=18,Ne=5,err;
  const double ei=8.14451245649177e-7,ef=322.538746367931;
  double logs,lsi,lsf,dls,*eos_t=NULL,z_table[5*18];
  eosp point;
  FILE *dadosout;
  
  lsi=-10.;//log(e2s_qgphr(ei,NULL));
  lsf= 6.5;//log(e2s_qgphr(ef,NULL));
    
  err= load_zoltan_table("zoltan.dat",z_table);
    
  dls=0.0005; // ds = 0.0005*hbarc;
  dadosout=fopen("zoltan_table2.eos","w");  
   
  for(logs = lsi ; logs < lsf ; logs += dls){
    point.s = exp(logs);
    err=zoltan_eos2(&point,z_table); 
    if(err!=0){
      continue;
    }
    fprintf(dadosout,"%.10lf %.10lf %.10lf %.10lf %.10lf\n",point.T,point.cs2,point.e,point.p,logs);
  }
  printf("end of loop\n");
  fclose(dadosout);
  
  if(eos_t!=NULL)
    free(eos_t);
  
  return 0;
}
