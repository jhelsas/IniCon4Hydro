/* Example: parse a simple configuration file */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ini.h"

typedef struct config
{
  int D,N_sph,N_species;
  double h,kh,ti,tf,dt,Tfo;
  double *xmin,*xmax,*dx;
  char *file,*hash,*folder,*integrator;
  char *kernel,*EoS, *setup, *deriv;
  char **sphp;
} config;

static int handler(void* user, const char* section, const char* name,
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
	
	for(l=0;l<(pcfg->D);l+=1){
      sprintf(keyname,"xmin[%d]",l);
      if(MATCH("acellerator-settings",keyname))
        pcfg->xmin[l]=atof(value);
        
      sprintf(keyname,"xmax[%d]",l);
      if(MATCH("acellerator-settings",keyname))
        pcfg->xmax[l]=atof(value);
        
      sprintf(keyname,"dx[%d]",l);
      if(MATCH("acellerator-settings",keyname))
        pcfg->dx[l]=atof(value);
    }
    
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

int print_cfg(config *cfg,FILE *cfgfile){
  int i,l;
  fprintf(cfgfile,"[simulation-parameters]\n");
  fprintf(cfgfile,"D = %d\n",cfg->D);
  fprintf(cfgfile,"N_sph = %d\n",cfg->D);
  fprintf(cfgfile,"N_species = %d\n",cfg->D);
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
  for(l=0;l<cfg->D;l+=1){
    fprintf(cfgfile,"xmin[%d] = %f\n",l,cfg->xmin[l]);
    fprintf(cfgfile,"xmax[%d] = %f\n",l,cfg->xmax[l]);
    fprintf(cfgfile,"dx[%d] = %f\n",l,cfg->dx[l]);
  }
  if(strcmp(cfg->file,"none")==0){
    fprintf(cfgfile,"\n[SPHparticle-data]\n");
    for(i=0;i< cfg->N_sph;i+=1)
      fprintf(cfgfile,"%d = %s\n",i,cfg->sphp[i]);
  }
  return 0;
}

int main(int argc, char* argv[])
{
  int i,l;
  FILE *outfile;
  config cfg;
  if(config_init(&cfg)!=0)
    return (-1);

  if (ini_parse("gubser.ini", handler, &cfg) < 0) {
    printf("Can't load 'test.ini'\n");
    return 1;
  }
  printf("Config loaded from gubser.ini: \n\nD=%d",cfg.D);
  printf(" N_sph = %d N_species=%d \nh=%lf kh=%lf\n",
           cfg.N_sph,cfg.N_species,cfg.h,cfg.kh);
  printf(" ti=%lf tf=%lf dt=%lf \n T_fo=%f\ndata_file=%s\n",
         cfg.ti,cfg.tf,cfg.dt,cfg.Tfo,cfg.file);
  printf("hash=%s\nfolder=%s\nkernel=%s\n",cfg.hash, cfg.folder,cfg.kernel);
  printf("EoS=%s\nsetup=%s\nderiv=%s\n",cfg.EoS, cfg.setup,cfg.deriv);
  printf("integrator = %s\n",cfg.integrator);
  printf("\n");
  for(l=0;l<(cfg.D);l+=1)
    printf("%d: %lf %lf %lf\n",l,cfg.xmin[l],cfg.xmax[l],cfg.dx[l]);       
  printf("\n");
  
  if( strcmp(cfg.file,"none") == 0)
    for(i=0;i < cfg.N_sph;i+=1)
      printf("%d: %s\n",i,cfg.sphp[i]);
      
  outfile=fopen("gubser_check.ini","w");
  if(print_cfg(&cfg,outfile) !=0)
    return 4;
  fclose(outfile);
      
  if(config_free(&cfg)!=0)
    return 3;
  return 0;
}
