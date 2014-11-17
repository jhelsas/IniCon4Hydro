#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include "inih_r29/ini.h"
#include "ini_rw.h"

int handler(void* user, const char* section, const char* name,
                   const char* value)
{
  int i,l;
  char keyname[100+1];
  config* pcfg = (config*)user;

  #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
  if (MATCH("simulation-parameters", "D")) 
    pcfg->D = atoi(value);
  else if (MATCH("simulation-parameters", "N_sph"))
    pcfg->N_sph = atoi(value);
  else if (MATCH("simulation-parameters", "N_species"))
    pcfg->N_species = atoi(value);
  else if (MATCH("simulation-parameters", "h"))
    pcfg->h = atof(value);
  else if (MATCH("simulation-parameters", "kh"))
    pcfg->kh = atof(value);
  else if (MATCH("simulation-parameters", "ti"))
    pcfg->ti = atof(value);
  else if (MATCH("simulation-parameters", "tf"))
    pcfg->tf = atof(value);
  else if (MATCH("simulation-parameters", "dt"))
    pcfg->dt = atof(value);
  else if (MATCH("simulation-parameters","T_fo"))
    pcfg->Tfo = atof(value);
  else if (MATCH("simulation-parameters", "data_file"))
    pcfg->file = strdup(value);
  else if (MATCH("simulation-parameters", "hash"))
    pcfg->hash = strdup(value);
  else if (MATCH("simulation-parameters","folder"))
    pcfg->folder = strdup(value);
  else if (MATCH("simulation-parameters","integrator"))
    pcfg->integrator = strdup(value);
  else if (MATCH("simulation-parameters","kernel"))
    pcfg->kernel = strdup(value);
  else if (MATCH("simulation-parameters","EoS"))
    pcfg->EoS = strdup(value);
  else if (MATCH("simulation-parameters","setup"))
    pcfg->setup = strdup(value);
  else if (MATCH("simulation-parameters","deriv"))
    pcfg->deriv = strdup(value);
  else if(pcfg->D>0){
	if(pcfg->xmin == NULL)
      pcfg->xmin = (double *)malloc((pcfg->D)*sizeof(double));
	
	if(pcfg->xmax == NULL)
      pcfg->xmax = (double *)malloc((pcfg->D)*sizeof(double));
	
	if(pcfg->dx == NULL)
      pcfg->dx  = (double *)malloc((pcfg->D)*sizeof(double));	    
	
	sprintf(keyname,"xmin");
    if(MATCH("acellerator-settings",keyname))
    for(l=0;l< pcfg->D;l+=1)
      sscanf(value,"%lf",&(pcfg->xmin[l]));  
	
	sprintf(keyname,"xmax");
    if(MATCH("acellerator-settings",keyname))
    for(l=0;l< pcfg->D;l+=1)
      sscanf(value,"%lf",&(pcfg->xmax[l]));
      
	sprintf(keyname,"dx");
    if(MATCH("acellerator-settings",keyname))
    for(l=0;l< pcfg->D;l+=1)
      sscanf(value,"%lf",&(pcfg->dx[l]));
    
    if(pcfg->N_sph > 0 && strcmp(pcfg->file,"none") == 0){
      if(pcfg->sphp==NULL)
        pcfg->sphp = (char**)malloc((pcfg->N_sph)*sizeof(char*));
      
      for(i=0;i < pcfg->N_sph;i+=1){
        sprintf(keyname,"%d",i);
        if(MATCH("SPHparticle-data",keyname))
          pcfg->sphp[i]=strdup(value);
      }
    }
  } else
    return 0;  /* unknown section/name, error */
  
  return 1;
}

int config_init(config*cfg){
	
  cfg->D=0; cfg->N_sph=0; cfg->sphp=NULL;
  cfg->Tfo=0.;cfg->ti=0.;cfg->tf=0.;cfg->dt=0.;
  cfg->xmax=NULL; cfg->xmin=NULL; cfg->dx=NULL; cfg->integrator=NULL;
  cfg->file=NULL; cfg->hash=NULL; cfg->folder=NULL; 
  cfg->kernel=NULL; cfg->EoS=NULL; cfg->setup=NULL; cfg->deriv=NULL;
  
  return 0;
}

int config_free(config *cfg){
  int i;
  if(strcmp(cfg->file,"none") == 0){
    for(i=0;i < cfg->N_sph;i+=1)
      free(cfg->sphp[i]);
    free(cfg->sphp);
  }
  if(cfg->file!=NULL)
    free(cfg->file); 
  if(cfg->hash!=NULL)
    free(cfg->hash);
  if(cfg->folder!=NULL) 
    free(cfg->folder);
  if(cfg->xmin!=NULL)
    free(cfg->xmin);
  if(cfg->xmax!=NULL) 
    free(cfg->xmax); 
  if(cfg->dx!=NULL)
    free(cfg->dx); 
  if(cfg->kernel!=NULL)
    free(cfg->kernel);
  if(cfg->EoS!=NULL)
    free(cfg->EoS);
  if(cfg->setup!=NULL)
    free(cfg->setup);
  if(cfg->deriv!=NULL)
    free(cfg->deriv);
  if(cfg->integrator!=NULL)
    free(cfg->integrator);
  
  return 0;
}

int fprint_cfg(config *cfg,FILE *cfgfile){
  int i,l;
  // write a comment preamble for users
  fprintf(cfgfile,"[simulation-parameters]\n");
  fprintf(cfgfile,"D = %d\n",cfg->D);
  fprintf(cfgfile,"N_sph = %d\n",cfg->N_sph);
  fprintf(cfgfile,"N_species = %d\n",cfg->N_species);
  fprintf(cfgfile,"kernel = %s\n",cfg->kernel);
  fprintf(cfgfile,"integrator = %s\n",cfg->integrator);
  fprintf(cfgfile,"EoS = %s\n",cfg->EoS);
  fprintf(cfgfile,"setup = %s\n",cfg->setup);
  fprintf(cfgfile,"deriv = %s\n",cfg->deriv);
  fprintf(cfgfile,"h = %f\n",cfg->h);
  fprintf(cfgfile,"kh = %f\n",cfg->kh);
  fprintf(cfgfile,"ti = %f\n",cfg->ti);
  fprintf(cfgfile,"tf = %f\n",cfg->tf);
  fprintf(cfgfile,"dt = %f\n",cfg->dt);
  fprintf(cfgfile,"T_fo = %f\n",cfg->Tfo);
  fprintf(cfgfile,"data_file = %s\n",cfg->file);
  fprintf(cfgfile,"hash = %s\n",cfg->hash);
  fprintf(cfgfile,"folder = %s\n",cfg->folder);
  
  fprintf(cfgfile,"\n[accelerator-settings]\n");
  if(cfg->xmin!=NULL){
    fprintf(cfgfile,"xmin = ");
    for(l=0;l< cfg->D;l+=1)
      fprintf(cfgfile,"%lf ",cfg->xmin[l]);}
  
  if(cfg->xmax!=NULL){
    fprintf(cfgfile,"\nxmax = ");
    for(l=0;l< cfg->D;l+=1)
      fprintf(cfgfile,"%lf ",cfg->xmax[l]);}
      
  if(cfg->dx!=NULL){  
    fprintf(cfgfile,"\ndx = ");
    for(l=0;l< cfg->D;l+=1)
      fprintf(cfgfile,"%lf ",cfg->dx[l]);}
    
  if(strcmp(cfg->file,"none")==0){
    fprintf(cfgfile,"\n\n[SPHparticle-data]\n");
    for(i=0;i< cfg->N_sph;i+=1)
      fprintf(cfgfile,"%d = %s\n",i,cfg->sphp[i]);
  }
  return 0;
}

int ini_simple_read(char* filename,int *Dout,int *Nout,
                    double **xout,double **uout,
                    double **Sout,config *cfg)
{
  int i,l,D,N;
  double ni,q,*x,*u,*S;
  
  if( config_init(cfg) != 0)
    return 1;
  
  if (ini_parse(filename, handler, cfg) < 0)
    return 2;
  
  D = cfg->D;
  N = cfg->N_sph;
  *Dout=D;
  *Nout=N;
  x = (double*)malloc(N*D*sizeof(double));
  u = (double*)malloc(N*D*sizeof(double));
  S = (double*)malloc(N*sizeof(double));
  
  if(strcmp(cfg->file,"none")==0){
    for(i=0;i<N;i+=1){
      sscanf(cfg->sphp[i],"%lf %lf %lf",&ni,&q,&(S[i]));
      for(l=0;l<D;l+=1)
        sscanf(cfg->sphp[i],"%lf ",&(x[i*D+l]));
      for(l=0;l<D;l+=1)
        sscanf(cfg->sphp[i],"%lf ",&(u[i*D+l]));
    }
  }
  else{
	FILE *infile;
	infile = fopen(cfg->file,"r");
    for(i=0;i<N;i+=1){
      fscanf(infile,"%lf %lf %lf",&ni,&q,&(S[i]));
      for(l=0;l<D;l+=1)
        fscanf(infile,"%lf ",&(x[i*D+l]));
      for(l=0;l<D;l+=1)
        fscanf(infile,"%lf ",&(u[i*D+l]));
    }
    fclose(infile);
  }
  *xout=x;
  *uout=u;
  *Sout=S;
    
  return 0;
}

int config_init_std(&cfg){
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
}
