#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "trial-functions.h"
#include "splitandfit.h"
#include <string.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;

domain :: domain(): xv(0){
  D=0;Nv=0;good=0;type=-1;
  S=0.;
}

domain :: domain(int dim, int Nvertex, int Type): xv(dim*Nvertex){
  D=dim; Nv=Nvertex; good=0; type=Type;
  S=0.;
}

int domain :: init(int dim, int Nvertex, int Type){
  D=dim; Nv=Nvertex; good=0; type=Type;
  xv.resize(dim*Nvertex);
  S=0.;
  return 0;
}

int domain_count(vector <domain> dom){
  int i=0;
  vector <domain> :: iterator el;
  for(el=dom.begin();el!=dom.end();el++)
    if(el->good==0)
      i+=1;
  return i;
}

int domain_check(vector <domain> dom){
  vector <domain> :: iterator el;
  for(el=dom.begin();el!=dom.end();el++)
    if(el->good!=0)
      return 1;  
  return 0;
}

int el_print(domain el){
  int j,l;
  for(j=0;j<el.Nv;j+=1){
    cout << "( ";
    for(l=0;l<el.D;l+=1)
      cout << el.xv[j*(el.D)+l] << " ";
    cout << ") ";
  }
  cout << endl;
  return 0;
}

int rollin(int *ind,int *indx,int D,int n){
  int i,index=0,tmp=1;
  if(D<=0 || n<=0)
    return 1;
  for(i=0;i<D;i+=1){
    index+=tmp*ind[i];tmp*=n;}
  *indx=index;
  return 0;
}

int rollout(int *ind,int *indx,int D,int n){
  int i,index=*indx;
  if(D<=0 || n<=0)
    return 1;
  for(i=0;i<D;i+=1){
    ind[i] = index%n;index=index/n;}
  return 0;
}

int cubic_split(int Nsplit,vector <domain> & dom, int n,int D){
  int i,j,k,l,err,indx,indi[D],ind[D],inda[D];
  double xv[(1<<D)*D];
  
  if(dom[n].type!=0)
    return -1;
  
  for(i=0;i<(1<<D);i+=1)
    for(l=0;l<D;l+=1)
      xv[i*D+l] = dom[n].xv[i*D+l];
    
  for(i=0;i<Nsplit;i+=1){
    domain tdel(D,1<<D,0); // vai tomar no cu, pq isso tem que estar aqui?
    tdel.good=0;
    err=rollout(ind,&i,D,2); if(err!=0) return err;
    
    for(l=0;l<D;l+=1){
      
      for(k=0;k<D;k+=1)
        inda[k] = ind[k];
      inda[l]=1-ind[l];
      
      err=rollin(inda,&indx,D,2); if(err!=0) return err;
      
      for(j=0;j<dom[n].Nv;j+=1){
        err=rollout(indi,&j,D,2); if(err!=0) return err;      
        if(indi[l]==ind[l])
          tdel.xv[D*j+l]=xv[D*i+l];
        else
          tdel.xv[D*j+l]=(xv[D*i+l]+xv[D*indx+l])/2.;
        
      }
    }
    tdel.S=0.0;
    dom.push_back(tdel);
  }
  
  return 0;
}

int init_cube(double xl[],double xu[],vector <domain> & dom,int D){
  int i,l,err,ind[D];
  domain mdel;
  err=mdel.init(D,1<<D,0); if(err!=0) return err;
  for(i=0;i<mdel.Nv;i+=1){
    err=rollout(ind,&i,D,2);
    for(l=0;l<D;l+=1){
      if(ind[l]==0)
        mdel.xv[i*D+l]=-1.;
      else
        mdel.xv[i*D+l]=1.;
    }
  }
  dom.push_back(mdel);
  return 0;
}

int clean_domain(vector <domain> &dom){
  unsigned int i;
  vector <domain> :: iterator el;  
  for(i=0;i<dom.size();i+=1)
    if(dom[i].good!=0){
      dom.erase(dom.begin()+i);
      i-=1;}
  if(domain_check(dom)!=0)
    cout<< "algum dominio esta ruim\n";
  
  return 0;
}

int print_sph(int D,const char *filename,vector <domain> dom){
  int l;
  ofstream sphfile;
  vector <domain> :: iterator el;
  if(D<=0 || filename==NULL || dom.size()==0)
    return 1;
  sphfile.open(filename);
  for(el=dom.begin();el!=dom.end();el++){
    for(l=0;l<D;l+=1) 
      sphfile << (el->xv[D*0+l]+el->xv[D*((1<<D)-1)+l])/2. << " ";
    for(l=0;l<D;l+=1) 
      sphfile << 0. << " ";
    sphfile << el->S <<"\n";
  }
  sphfile.close();
  return 0;
}

int domain_split(int D, int Nsplit,double cutoff,vector <domain>& dom, gsl_monte_function F){
  unsigned int id;
  int err,l;
  double xi[D],xf[D],S,erd;
  size_t calls = 500000;
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_monte_miser_state *s_m=gsl_monte_miser_alloc (D);
    
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  id=0;
  while(id<dom.size()){
    
    S=0.0;
    for(l=0;l<D;l+=1){
      xi[l]=dom[id].xv[D*0+l];
      xf[l]=dom[id].xv[D*((1<<D)-1)+l]; 
    }
    
    gsl_monte_miser_integrate(&F,xi,xf,D,calls,r,s_m,&S,&erd);
    
    dom[id].S=S;

    if(S > cutoff){
      dom[id].good=1;
      err=cubic_split(Nsplit,dom,id,D);if(err!=0){return err;}
    }
       
    id++;
  }
  
  gsl_monte_miser_free (s_m);
  gsl_rng_free (r);
  
  return 0;
}
