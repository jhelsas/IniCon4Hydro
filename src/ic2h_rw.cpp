#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include "inih_r29/ini.h"
#include "ini_rw.h"
#include "splitandfit.h"

using namespace std;

int ini_write(int D,const char *filename, char *hash, char *folder,
              vector <domain> &dom,
              int (*velocity)(double *,size_t,void *,double *),
              void *par, config *cfg)
{
  int err=0,j,l,it=0,N;
  double x[D],u[D];
  char part[200+1],buf[100+1];
  FILE *outfile;
  ofstream sphfile;
  vector <domain> :: iterator el;
  if(D<=0 || dom.size()==0)
    return 1;
  N=dom.size();
  cfg->D=D;
  cfg->N_sph = N;
  cfg->N_species = 0;
  cfg->kernel = strdup("b_spline");
  cfg->integrator = strdup("heun");
  cfg->EoS = strdup("qg");
  cfg->setup = strdup("ssph_2p1bi");
  cfg->deriv = strdup("Drv_2p1bi");
  cfg->hash = strdup(hash);
  cfg->folder = strdup(folder);
  cfg->h = 0.1;
  cfg->kh = 0.1;
  cfg->ti = 1.0;
  cfg->tf = 8.0;
  cfg->dt = 0.025;
  cfg->Tfo = 0.145;
  
  cfg->xmin=(double*)malloc(D*sizeof(double));
  cfg->xmax=(double*)malloc(D*sizeof(double));
  cfg->dx=(double*)malloc(D*sizeof(double));
  
  for(l=0;l<D;l+=1){
    cfg->xmin[l]=-20.; cfg->xmax[l]=20.; cfg->dx[l]=0.2;
  }  
  
  if(strcmp(filename,"none") == 0){
    cfg->file = strdup("none");
    cfg->sphp = (char **)malloc(N*sizeof(char*));
    for(el=dom.begin();el!=dom.end();el++){
      sprintf(part,"1.0 1.0 %lf ",el->S);
      if(el->type==0){
        for(l=0;l<D;l+=1){ 
          x[l] = (el->xv[D*0+l]+el->xv[D*((1<<D)-1)+l])/2.;
		  sprintf(buf,"%lf ",x[l]);
          strcat(part,buf);
        }
      }
      else{
        for(l=0;l<D;l+=1){
          x[l]=0.;
          for(j=0;j<el->Nv;j+=1)
            x[l]+= el->xv[D*j+l];
          x[l]=x[l]/(el->Nv);
          sprintf(buf,"%lf ",x[l]);
          strcat(part,buf);
        }
      }
      err=velocity(x,D,par,u);
      for(l=0;l<D;l+=1){
		sprintf(buf,"%lf ",u[l]);
        strcat(part,buf);
	  }
      cfg->sphp[it] = strdup(part);
      it++;
    }
  }
  else{
	cfg->file = strdup(filename);
    sphfile.open(filename);
    for(el=dom.begin();el!=dom.end();el++){
      sphfile << 1.0 << " " << 1.0 << " " << el->S << " ";
      if(el->type==0){
        for(l=0;l<D;l+=1){ 
          x[l] = (el->xv[D*0+l]+el->xv[D*((1<<D)-1)+l])/2.;
          sphfile << x[l] << " ";
        }
      }
      else{
        for(l=0;l<D;l+=1){
          x[l]=0.;
          for(j=0;j<el->Nv;j+=1)
            x[l]+= el->xv[D*j+l];
          x[l]=x[l]/(el->Nv);
          sphfile << x[l] << " ";
        }
      }
      err=velocity(x,D,par,u);
      for(l=0;l<D;l+=1){
        sphfile << u[l] << " ";
      }
      sphfile << "\n";
      it++;
    }
    sphfile.close();
  }
  
  if(strcmp(hash,"") !=0)
    sprintf(part,"%s.cfg",hash);
  else
    sprintf(part,"hydro.cfg");
    
  outfile = fopen(part,"w");
  err=fprint_cfg(cfg,outfile);
  fclose(outfile);
  
  return 0;
}
