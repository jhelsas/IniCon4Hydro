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
  
  dwds = (par->cs2+1.)*(par->T);
  par->hsh=(4.0/3.0)*(par->p); /* verificar isso */
    
  return 0;
}


gsl_interp_accel *acc; // global variables
gsl_spline *splT, *splp ,*sple ,*splcs2 ;

int eos_zoltan(eosp *par, void *params){
  static int call_count = 0;
  const int Ne=5,Np=18;
  const double hbarc=0.1973269718;
  const double ls[18] = {-1.769049, -1.077466, -0.395522, 0.142758, 0.549412, 0.824491, 1.108151, 1.425858, 1.732869, 2.388880, 2.923345, 3.268573, 3.909100, 4.582473, 5.569299, 6.140774, 7.031874, 7.722272 };
  const double T[18]  = {0.100000, 0.115000, 0.129000, 0.139000, 0.147000, 0.152000, 0.158000, 0.166000, 0.175000, 0.200000, 0.228000, 0.250000, 0.299000, 0.366000, 0.500000, 0.600000, 0.800000, 1.000000 };
  const double e[18]  = {0.014186, 0.032551, 0.073524, 0.137981, 0.221213, 0.302902, 0.419333, 0.602841, 0.858120, 1.844991, 3.499477, 5.333056, 11.848111, 27.884914, 100.540059, 212.359345, 686.086922, 1707.554136 };
  const double p[18]  =  {0.002863, 0.006601, 0.013335, 0.022349, 0.033425, 0.043768, 0.059210, 0.087956, 0.131831, 0.335264, 0.742100, 1.235398, 3.058248, 7.893719, 30.585002, 66.288501, 219.633110, 550.530030 };
  const double cs2[18]=  {0.190000, 0.180000, 0.140000, 0.130000, 0.120000, 0.120000, 0.140000, 0.160000, 0.180000, 0.220000, 0.260000, 0.270000, 0.290000, 0.320000, 0.320000, 0.320000, 0.320000, 0.320000 };
  double logs;
  
  if(call_count==0){
    acc=gsl_interp_accel_alloc ();
    splT = gsl_spline_alloc (gsl_interp_cspline, 18);
    splp = gsl_spline_alloc (gsl_interp_cspline, 18);
    sple = gsl_spline_alloc (gsl_interp_cspline, 18);
    splcs2 = gsl_spline_alloc (gsl_interp_cspline, 18);
    gsl_spline_init (splT, ls, T, Np);
    gsl_spline_init (splp, ls, p, Np);
    gsl_spline_init (sple, ls, e, Np);
    gsl_spline_init (splcs2, ls, cs2, Np);  
  }
  call_count+=1;
  
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

double e2s_zoltanRF(double s,void *params){
  int err;
  double *e0=(double*)params;
  eosp par;
  par.s=s;
  err = eos_zoltan(&par,NULL);
  
  printf("%lf %lf %lf %lf %lf\n",par.s,par.T,par.e,par.p, par.cs2);
  printf("s=%lf e =%lf e0 =%lf \n",s,par.e,*e0);
  scanf("%d",&err);
  
  return (par.e - *e0);
}

const gsl_root_fsolver_type *rsT;
gsl_root_fsolver *solver=NULL;
    
double e2s_zoltan(double epsilon, void *p){
  static int call_count = 0;
  int status, iter, max_iter;
  const double lsi = -2.3, lsf=6.5, ei=0.014186, ef=1707.554136;
  double x_lo,x_hi;
  gsl_function F;
  double s;
    
  if(call_count==0){
    rsT = gsl_root_fsolver_brent;
    solver = gsl_root_fsolver_alloc (rsT);
  }
  call_count +=1;
  
  F.function = &e2s_zoltanRF;
  F.params = &epsilon;
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
  
  return s;
}

int main(){
  int i,Np=18,Ne=5,err;
  const double ei=0.1,ef=322.538746367931;
  double e,de=0.005,logs,lsi,lsf,dls,*eos_t=NULL,z_table[5*18];
  eosp point;
  FILE *dadosout;
  
  lsi=-2.3;//log(e2s_qgphr(ei,NULL));
  lsf= 6.5;//log(e2s_qgphr(ef,NULL));
    
  err= load_zoltan_table("zoltan.dat",z_table);
  if(err!=0){printf("zoltan table not loaded\n");return err;}
  
  dls=0.0005; // ds = 0.0005*hbarc;
  dadosout=fopen("zoltan_table.eos","w");  
   
  for(logs = lsi ; logs < lsf ; logs += dls){
    point.s = exp(logs);
    err=eos_zoltan(&point,z_table); 
    //err=eos_zoltan2(&point,z_table); 
    //err=eos_zoltan3(&point,z_table); 
    if(err!=0)
      continue;
    fprintf(dadosout,"%.10lf %.10lf %.10lf %.10lf %.10lf\n",point.T,point.cs2,point.e,point.p,logs);
  }
  printf("end of loop\n");
  fclose(dadosout);
  /*
  dadosout = fopen("zoltan_e2s_comparison.dat","w");
  for(e=ei;e<=ef;e+=de){
    logs = e2s_zoltan(e,NULL);
    point.s = logs;
    err=eos_zoltan(&point,NULL); 
    if(err!=0)
      continue;
    fprintf(dadosout,"%lf %lf %lf\n",logs,e,point.e);
  }
  fclose(dadosout);
  */
  if(eos_t!=NULL)
    free(eos_t);
    /*
  gsl_spline_free(splT); 
  gsl_spline_free(sple); 
  gsl_spline_free(splp); 
  gsl_spline_free(splcs2);
  gsl_interp_accel_free(acc); */
  
  //if(solver!=NULL)
  //  gsl_root_fsolver_free (solver);
  
  return 0;
}
