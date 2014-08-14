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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;

/*
 * Arquivo de implementação das rotinas de divisão e fit
 * 
 * Leia splitandfit.h para quais funções usar e quais funções 
 * são internas da biblioteca
 * 
 */

/*
 * Public Functions
 * 
 * domain():
 *   Default Constructor
 * 
 * domain(int,int,int):
 *   Constructor
 * 
 * init(int,int,int):
 *  Initializer Fuction
 *  Usage: Starting From the Default Constructor,
 *    Input the dimension, Number of vertex and type of domain 
 */

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

/*
 * Public Functions
 * 
 * domain_count: 
 *   Output the number of domains on dom
 * 
 * domain_check:
 *   check if all the domains are good, or if there are spurious domains
 * 
 * el_print:
 *   Print one domain coordinates 
 * 
 * print_domain_n:
 *   Debug function - print the n-th element of dom's vertex coordinates
 * 
 * print_domain:
 *   print the coordinates of all vertexes of the domain list
 * 
 * check_inside:
 *   check if a point is inside a given domain
 *   to use with functions as input for the code
 */
 
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

int print_domain_n(int D,vector <domain> dom,int n){
  int i,l;
 
  cout << "domain " << n << ":\n";
  for(i=0;i<dom[n].Nv;i+=1){
    cout << "  vertex " << i << ": ";
    for(l=0;l<D;l+=1)
      cout << dom[n].xv[i*D+l] << " ";
    cout << "\n";
  }
  return 0;
}

int print_domain(int D,vector <domain> dom){
  int n,i,l;
  for(n=0;n<dom.size();n+=1){
    cout << "domain " << n << ":\n";
    for(i=0;i<dom[n].Nv;i+=1){
      cout << "  vertex " << i << ": ";
      for(l=0;l<D;l+=1)
        cout << dom[n].xv[i*D+l] << " ";
      cout << "\n";
    }
    cout << "aspect ratio: " << aspect_ratio(dom[n]) << "\n";
  }
  cout << "\n";
  return 0;
}

int check_inside(double *x,size_t D,domain *mdel){
  int i,l,err;
  double lmb[D],min=100000.,max=-1000000.;
  
  if(mdel->type == 0){
    return 0;
  }
  else if(mdel->type == 1){
    err=bc_coord(D,x,lmb,*mdel);if(err!=0){printf("problemas gsl_wrapper\n");return 0;}
    err=bc_check(D,lmb);
    if(err==0)
      return 0;
    else
      return 1;
  }
  else
    return 1;
  return 1;
}

/*
 * Public Function
 * 
 * init_cube:
 *   Initalize the domain with a (d-)cube domain aligned with 
 *   artesian axis, the lower and upper limits of the cube are given by
 *   xl[] and xu[] 
 * 
 *   e.g.:
 *   2-cube = square (regular) or paralelogram (irregular)
 *   3-cube = cube (regular) or paralelepiped (irregular)
 * 
 * init_triangle:
 *   
 *   Initalize the domain with a (d-)tiangles/simplex form the number of
 *   triangles and from the list of the triangles positions
 * 
 *   e.g.:
 *   2-triangle/simplex = triangles (regular or irregular)
 *   3-triangle/simples = tetrahedron (regular or irregular)
 */

int init_cube(double xl[],double xu[],vector <domain> & dom,int D){
  int i,l,err,ind[D];
  domain mdel;
  if(dom.size()>0)
    return 1;
  err=mdel.init(D,1<<D,0); if(err!=0) return err;
  for(i=0;i<mdel.Nv;i+=1){
    err=rollout(ind,&i,D,2);
    for(l=0;l<D;l+=1){
      if(ind[l]==0)
        mdel.xv[i*D+l]=xl[l]; 
      else
        mdel.xv[i*D+l]=xu[l];
    }
  }
  dom.push_back(mdel);
  return 0;
}

int init_triangle(int D,double Ntri,double *xv,vector <domain> & dom){
  int i,n,l,err;
  if(dom.size()>0)
    return 1;
  
  for(n=0;n<Ntri;n+=1){
    domain mdel;
    err=mdel.init(D,D+1,1); if(err!=0) return err;
    for(i=0;i<mdel.Nv;i+=1)
      for(l=0;l<D;l+=1)
        mdel.xv[i*D+l] = xv[n*D*(D+1)+i*D+l];
    dom.push_back(mdel);
  }

  return 0;
}

/*
 * Public Function
 * 
 * domain_split:
 *   Main function of the code
 *   
 *   Takes a domain list dom, and splits til the weight within it 
 *   non-breaked domain is less than a given cutoff
 * 
 *   The weight is calculated as the integral of the function F
 *   on the associated domain. The breaking is done automatically, and the
 *   function returns either sucess or fail.
 * 
 *   The result is stored on the input domain list dom. It stills contains
 *   spurious (breaked), thus, it's necessary to call clean_domain before
 *   turning it into SPH particles
 */
 
int domain_split(int D,double cutoff,vector <domain>& dom, gsl_monte_function F){
  unsigned int id;
  int i,err,l;
  double xl[D],xu[D],S,erd;
  size_t calls = 500000;
  const gsl_rng_type *T;
  wparams *lpar;
  gsl_rng *r;
  gsl_monte_miser_state *s_m=gsl_monte_miser_alloc (D);
    
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  lpar=(wparams*)(F.params);
    
  id=0;
  while(id<dom.size()){
    
    lpar->mdel=&(dom[id]);
    
    S=0.0;
    if(dom[id].type==0){
      for(l=0;l<D;l+=1){
        xl[l]=dom[id].xv[D*0+l];
        xu[l]=dom[id].xv[D*((1<<D)-1)+l]; 
      }
    }
    else{
      for(l=0;l<D;l+=1){
        xl[l]=100000.;
        xu[l]=-10000.;
        for(i=0;i<(dom[id].Nv);i+=1){
          if(dom[id].xv[D*i+l] < xl[l])
            xl[l] = dom[id].xv[D*i+l];
          if(dom[id].xv[D*i+l] > xu[l])
            xu[l] = dom[id].xv[D*i+l];
        }
      }
    }
    
    gsl_monte_miser_integrate(&F,xl,xu,D,calls,r,s_m,&S,&erd);
    
    dom[id].S=S;
    
    if(S > cutoff){
      dom[id].good=1;
      if(dom[id].type==0){
        err=cubic_split(D,dom,id);if(err!=0){return err;}}
      else if(dom[id].type==1){
        if(D==2){
          err=triangle_midpoint_split(D,dom,id);if(err!=0){ cout << "in: " << err << endl;return err;}}
        else{
          err=bc_simplex_split(D,dom,id);if(err!=0) return err;}
      }
    }
    id++;
  }
  
  gsl_monte_miser_free (s_m);
  gsl_rng_free (r);
  
  return 0;
}

/*
 * Public Functions
 * 
 * clean_domain:
 *   Use after domain_split
 *   Clean Spurious (already split) Domains
 * 
 * print_sph:
 *   Use after clean_domain
 *   Converts the domains weights and positions to produce 
 *   SPH particles
 */

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

int print_sph(int D,const char *filename,vector <domain> &dom){
  int j,l,it=0;
  ofstream sphfile;
  vector <domain> :: iterator el;
  if(D<=0 || filename==NULL || dom.size()==0)
    return 1;
  sphfile.open(filename);
  for(el=dom.begin();el!=dom.end();el++){
    if(el->type==0){
      for(l=0;l<D;l+=1) 
        sphfile << (el->xv[D*0+l]+el->xv[D*((1<<D)-1)+l])/2. << " ";
      for(l=0;l<D;l+=1) 
        sphfile << 0. << " ";
    }
    else{
      double x[l];
      for(l=0;l<D;l+=1){
        x[l]=0.;
        for(j=0;j<el->Nv;j+=1)
          x[l]+= el->xv[D*j+l];
        x[l]=x[l]/(el->Nv);
        sphfile << x[l] << " ";
      }
      for(l=0;l<D;l+=1)
        sphfile << 0. << " ";
    }
    sphfile << el->S <<"\n";
    it++;
  }
  sphfile.close();
  return 0;
}

/* 
 * Internal Function
 * 
 * rollin , rollout , aspect_ratio
 */

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

double aspect_ratio(domain mdel){
  int i,j,l;
  double dist,min=100000000.,max=-1000000000.;
  if(mdel.type!=1)
    return -1.;
  for(i=0;i<mdel.Nv;i+=1){
    for(j=0;j<mdel.Nv;j+=1){
      if(i==j)
        continue;
      dist=0.;
      for(l=0;l<mdel.D;l+=1)
        dist+= (mdel.xv[i*(mdel.D)+l]-mdel.xv[j*(mdel.D)+l])*(mdel.xv[i*(mdel.D)+l]-mdel.xv[j*(mdel.D)+l]);
      dist=sqrt(dist);
      if(dist>max)
        max=dist;
      if(dist<min)
        min=dist;
    }
  }
  return max/min;
}

/*
 * Internal Function
 * 
 * Spliting routine for (d-dimensional) cubic domains
 * Split 1 (d-)cube into 2^d smaller (d-)cubes
 * 
 * e.g.:
 * 2-cube = square (regular) or paralelogram (irregular)
 * 3-cube = cube (regular) or paralelepiped (irregular)
 */

int cubic_split(int D,vector <domain> & dom, int n){
  int i,j,k,l,err,indx,indi[D],ind[D],inda[D];
  double xv[(1<<D)*D];
  
  if(dom[n].type!=0)
    return -1;
  
  for(i=0;i<(1<<D);i+=1)
    for(l=0;l<D;l+=1)
      xv[i*D+l] = dom[n].xv[i*D+l];
    
  for(i=0;i<(1<<D);i+=1){
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

/*
 * Internal function
 * 
 * Spliting routine for (d-dimensional) simplex into d+1 simplexes 
 * using the baricenter as a vertex. 
 * Note: Not optimal split
 * 
 * e.g.:
 * 2-triangle/simplex = triangles (regular or irregular)
 * 3-triangle/simples = tetrahedron (regular or irregular)
 */

int bc_simplex_split(int D,vector <domain> & dom,int n){
  int i,j,l,ind; 
  double bc[D];
  
  if(dom[n].type!=1)
    return -1;
  for(l=0;l<D;l+=1){
    bc[l]=0.;
    for(i=0;i<dom[n].Nv;i+=1)
      bc[l]+=dom[n].xv[D*i+l];
    bc[l]*= 1./((double)dom[n].Nv);
  }
  
  for(i=0;i<dom[n].Nv;i+=1){
    domain tdel(D,D+1,1);
    tdel.good=0;
    for(j=0;j<=D;j+=1){
      for(l=0;l<D;l+=1){
        if(j==i)
          tdel.xv[D*j+l]=bc[l];
        else
          tdel.xv[D*j+l]=dom[n].xv[D*j+l];
      }
    }
    tdel.S=0.0;
    dom.push_back(tdel);
  }
  return 0;
}

/*
 * Internal functions
 * 
 * bc_coord:
 * Calculate the baricentric coordinates *lmb of a given point *r 
 * to a given (d-) triangular domain
 *  
 * bc_check
 * Check if a point is inside in a given (d-)triangular domain
 * using it's baricentric coordinates.
 * 
 * bc_test:
 * print the baricentric coordinates of *r and if they are inside
 * the domain or not
 */

int bc_coord(int D,double *r,double *lmb,domain mdel){
  int l,n,s;
  double T[D*D],B[D];
  if(mdel.type!=1)
    return -1;
  
  for(l=0;l<D;l+=1){
    B[l] = r[l] - mdel.xv[D*D+l];
    for(n=0;n<D;n+=1)
      T[l*D+n] = mdel.xv[n*D+l] - mdel.xv[D*D+l];
  }
  
  gsl_matrix_view m = gsl_matrix_view_array (T, D, D);
  gsl_vector_view b = gsl_vector_view_array (B, D);
  gsl_vector *x = gsl_vector_alloc (D);  
  gsl_permutation * p = gsl_permutation_alloc (D);
  gsl_linalg_LU_decomp (&m.matrix, p, &s);
  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
  for(l=0;l<D;l+=1)
    lmb[l] = gsl_vector_get(x,l);
  gsl_permutation_free (p);
  gsl_vector_free (x);
  
  return 0;
}
 
int bc_check(int D,double *lmb){
  int l;
  double Lmb;
  if(D<=0)
    return -1;
  Lmb=0.;
  for(l=0;l<D;l+=1){
    Lmb+=lmb[l];
    if(lmb[l]<0 || lmb[l]>1)
      return 1;
  }
  if(Lmb <0 || Lmb > 1)
    return 1;
  return 0;
}

int bc_test(int D,double *r,domain mdel){
  int err,l;
  double lmb[D],Lmb;
  err=bc_coord(D,r,lmb,mdel);if(err!=0) return err;
  cout << "bc coordinates: " ;
  Lmb=0.;
  for(l=0;l<D;l+=1){
    cout << lmb[l] << " " ;Lmb+= lmb[l];}
  cout << Lmb << "\n ";
    
  err= bc_check(D,lmb);
  if(err==0)
    cout << "ponto no interior\n";
  else
    cout << "ponto no exterior\n";
  return 0;
}
/*
 * Internal function
 * 
 * Spliting routine for (2-dimensional) triangles into 4 triangles 
 * midpoints of each side as vertexes.
 * 
 */
 
int triangle_midpoint_split(int D,vector <domain> & dom,int n){
  int i,l,ind1,ind2;
  double midp[(D+1)*D];
  
  if(dom[n].type!=1 || D!=2)
    return -1;
  
  for(i=0;i<dom[n].Nv;i+=1){
    ind1=(i+1)%(dom[n].Nv);
    ind2=(i+2)%(dom[n].Nv);
    for(l=0;l<D;l+=1)
      midp[i*D+l] = (dom[n].xv[ind1*D+l]+dom[n].xv[ind2*D+l])/2.;
  }
  
  for(i=0;i<dom[n].Nv;i+=1){
    domain tdel(D,D+1,1);
    tdel.good=0;
    
    ind1=(i+1)%(dom[n].Nv);
    ind2=(i+2)%(dom[n].Nv);
    
    for(l=0;l<D;l+=1){
      tdel.xv[D*i+l]=dom[n].xv[D*i+l];
      tdel.xv[D*ind1+l]=midp[ind1*D+l];
      tdel.xv[D*ind2+l]=midp[ind2*D+l];
    }
    tdel.S=0.0;
    dom.push_back(tdel);
  }
  
  {
    domain tdel(D,D+1,1);
    tdel.good=0;
    
    for(i=0;i<dom[n].Nv;i+=1)
      for(l=0;l<D;l+=1)
        tdel.xv[D*i+l] = midp[i*D+l];
        
    tdel.S=0.0;
    dom.push_back(tdel);
  }
  
  return 0;
}

/*
 * Internal function
 * 
 * Spliting routine for (3-dimensional) tetrahedrons into 8 triangles 
 * midpoints of each side as vertexes.
 * 
 */
 
/*
int tetrahedron_midpoint_split(int D,vector <domain> & dom,int n){
  int i,j,l,ind[n];
  double midp[(D+1)*D*D];
  
  if(dom[n].type!=1 || D!=3)
    return -1;
    
  for(i=0;i<dom[n].Nv;i+=1){
    for(j=i;j<dom[n].Nv;j+=1)
      ind[i]=(i+j)%(dom[n].Nv);
    
    for(j=i;j<dom[n].){
      for(l=0;l<D;l+=1)
        midp[((D+1)*i+j)*D+l] = 
    }
  }
  
  return 0;
}*/


/*
 * Public Functions
 * 
 * unit_hexagon and unit2_hexagon
 *   Utility Functions. Create 6 triangles that form an hexagon
 *   to use as input for init_triangle with Ntri=6
 * 
 *   unit_hexagon creates an hexagon with side = 1.
 *   unit2_hexagon creates an hexagon with side = 2.
 */

int unit_hexagon(int Ntri,int D,double *xv){
  if(Ntri!=6 || D!=2)
    return -1;
    
  xv[0*(D+1)*D+0*D+0]=1.; xv[0*(D+1)*D+0*D+1]=0.;
  xv[0*(D+1)*D+1*D+0]=.5; xv[0*(D+1)*D+1*D+1]=sqrt(3.)/2.;
  xv[0*(D+1)*D+2*D+0]=0.; xv[0*(D+1)*D+2*D+1]=0.;
  
  xv[1*(D+1)*D+0*D+0]=.5; xv[1*(D+1)*D+0*D+1]=sqrt(3.)/2.;
  xv[1*(D+1)*D+1*D+0]=-.5;xv[1*(D+1)*D+1*D+1]=sqrt(3.)/2.;
  xv[1*(D+1)*D+2*D+0]=0.; xv[1*(D+1)*D+2*D+1]=0.;
  
  xv[2*(D+1)*D+0*D+0]=0.; xv[2*(D+1)*D+0*D+1]=0.;
  xv[2*(D+1)*D+1*D+0]=-.5;xv[2*(D+1)*D+1*D+1]=sqrt(3.)/2.;
  xv[2*(D+1)*D+2*D+0]=-1.;xv[2*(D+1)*D+2*D+1]=0.;
  
  xv[3*(D+1)*D+0*D+0]=-1.;xv[3*(D+1)*D+0*D+1]=0.;
  xv[3*(D+1)*D+1*D+0]=-.5;xv[3*(D+1)*D+1*D+1]=-sqrt(3.)/2.;
  xv[3*(D+1)*D+2*D+0]=0.; xv[3*(D+1)*D+2*D+1]=0.;
  
  xv[4*(D+1)*D+0*D+0]=-.5;xv[4*(D+1)*D+0*D+1]=-sqrt(3.)/2.;
  xv[4*(D+1)*D+1*D+0]=.5; xv[4*(D+1)*D+1*D+1]=-sqrt(3.)/2.;
  xv[4*(D+1)*D+2*D+0]=0.; xv[4*(D+1)*D+2*D+1]=0.;
  
  xv[5*(D+1)*D+0*D+0]=.5;xv[5*(D+1)*D+0*D+1]=-sqrt(3.)/2.;
  xv[5*(D+1)*D+1*D+0]=1.;xv[5*(D+1)*D+1*D+1]=0.;
  xv[5*(D+1)*D+2*D+0]=0.;xv[5*(D+1)*D+2*D+1]=0.;
  
  return 0;
}
 
int unit2_hexagon(int Ntri,int D,double *xv){
  if(Ntri!=6 || D!=2)
    return -1;
    
  xv[0*(D+1)*D+0*D+0]=2.; xv[0*(D+1)*D+0*D+1]=0.;
  xv[0*(D+1)*D+1*D+0]=1.; xv[0*(D+1)*D+1*D+1]=sqrt(3.);
  xv[0*(D+1)*D+2*D+0]=0.; xv[0*(D+1)*D+2*D+1]=0.;
  
  xv[1*(D+1)*D+0*D+0]=1.; xv[1*(D+1)*D+0*D+1]=sqrt(3.);
  xv[1*(D+1)*D+1*D+0]=-1.;xv[1*(D+1)*D+1*D+1]=sqrt(3.);
  xv[1*(D+1)*D+2*D+0]=0.; xv[1*(D+1)*D+2*D+1]=0.;
  
  xv[2*(D+1)*D+0*D+0]=0.; xv[2*(D+1)*D+0*D+1]=0.;
  xv[2*(D+1)*D+1*D+0]=-1.;xv[2*(D+1)*D+1*D+1]=sqrt(3.);
  xv[2*(D+1)*D+2*D+0]=-2.;xv[2*(D+1)*D+2*D+1]=0.;
  
  xv[3*(D+1)*D+0*D+0]=-2.;xv[3*(D+1)*D+0*D+1]=0.;
  xv[3*(D+1)*D+1*D+0]=-1.;xv[3*(D+1)*D+1*D+1]=-sqrt(3.);
  xv[3*(D+1)*D+2*D+0]=0.; xv[3*(D+1)*D+2*D+1]=0.;
  
  xv[4*(D+1)*D+0*D+0]=-1.;xv[4*(D+1)*D+0*D+1]=-sqrt(3.);
  xv[4*(D+1)*D+1*D+0]=1.; xv[4*(D+1)*D+1*D+1]=-sqrt(3.);
  xv[4*(D+1)*D+2*D+0]=0.; xv[4*(D+1)*D+2*D+1]=0.;
  
  xv[5*(D+1)*D+0*D+0]=1.;xv[5*(D+1)*D+0*D+1]=-sqrt(3.);
  xv[5*(D+1)*D+1*D+0]=2.;xv[5*(D+1)*D+1*D+1]=0.;
  xv[5*(D+1)*D+2*D+0]=0.;xv[5*(D+1)*D+2*D+1]=0.;
  
  return 0;
}
