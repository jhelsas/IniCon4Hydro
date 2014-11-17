/* Example: parse a simple configuration file */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "inih_r29/ini.h"
#include "ini_rw.h"

using namespace std;

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
