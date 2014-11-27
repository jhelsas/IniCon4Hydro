/*
 * Ideal massless EoS
 * 
 * e = C s^{4/3}
 */
 
#include <cmath>
#include "conv_funct.h"

/*
 * Tem que rever isso depois
 */

double e2s_pion(double epsilon,void *p){
  const double C_pi = 0.134828384; /* 3*(hbarc)*((45÷(3x128×π^2))^(1/3)) GeV fm */
  return pow(epsilon/(3.0*C_pi),3./4.);
}

/*
 * 3 quark flavor massless qgp
 */

double e2s_qg(double epsilon,void *p){
  //const double C_qg = 0.01948439; /* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */  
  const double C_qg = 0.0179276286445;
  return pow(epsilon/(3.0*C_qg),3./4.);
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


gsl_interp_accel *acc; // global variables
gsl_spline *splT, *splp ,*sple ,*splcs2 ;

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
    
double e2s_zoltan(double e, void *p){
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
    
  if(e < ei)
    return (e/ei)*0.170495053; // 0.170495053 * (e/ei)
    
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

double e2s_table(double epsilon,void *p){
  int i,Ne=5,Np,c;
  double *eos_t,logs;
  eos_t=(double*)p;
  
  if(epsilon > eos_t[Ne*(Np-1)+3] || epsilon < eos_t[Ne*0+3])
    return -1.;
  
  c=0;
  i=Np/2;
  while(c==0){
    if( epsilon < eos_t[(i+1)*Ne+3] )
      i = i/2;
    else if( epsilon > eos_t[(i+1)*+3] )
      i = (i+Ne)/2;
    else
      c=1;
  }
  
  logs = (epsilon - eos_t[i*Ne+3])/(eos_t[(i+1)*Ne+3]-eos_t[i*Ne+3]);
  logs = logs*(eos_t[(i+1)*Ne+0]-eos_t[i*Ne]) + eos_t[i*Ne];
    
  return exp(logs);
}

double T2s_table(double T,void *p){
  int i,Ne=5,Np,c;
  double *eos_t;
  eos_t=(double*)p;
  
  if(T > eos_t[Ne*(Np-1)+1] || T < eos_t[Ne*0+1])
    return -1.;
  
  c=0;
  i=Np/2;
  while(c==0){
    if( T < eos_t[(i+1)*Ne+3] )
      i = i/2;
    else if( T > eos_t[(i+1)*+3] )
      i = (i+Ne)/2;
    else
      c=1;
  }
  
  T = (T - eos_t[i*Ne+1])/(eos_t[(i+1)*Ne+1]-eos_t[i*Ne+1]);
  T = T*(eos_t[(i+1)*Ne]-eos_t[i*Ne]) + eos_t[i*Ne];
    
  return T;
}
