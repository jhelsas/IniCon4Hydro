typedef struct config
{
  int D,N_sph,N_species;
  double h,kh,ti,tf,dt,Tfo;
  double *xmin,*xmax,*dx;
  char *file,*hash,*folder,*integrator;
  char *kernel,*EoS, *setup, *deriv;
  char **sphp;
} config;

int handler(void* user, const char* section, const char* name,
                   const char* value);

int config_init(config*cfg);

int config_free(config *cfg);

int fprint_cfg(config *cfg,FILE *cfgfile);
