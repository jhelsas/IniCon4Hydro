/*
 * Ideal massless EoS
 * 
 * e = C s^{4/3}
 */
 
#include <cmath>

double e2s_pion(double epsilon,void *p){
  const double C_pi = 0.134828384; /* 3*(hbarc)*((45÷(3x128×π^2))^(1/3)) GeV fm */
  return pow(epsilon/C_pi,3./4.);
}

double e2s_qg(double epsilon,void *p){
  const double C_qg = 0.058356312;/* 3*(hbarc)*((45÷(37x128×π^2))^(1/3)) GeV fm */ 
  return pow(epsilon/C_qg,3./4.);
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
