#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int load_eos_table(char *filename,double *eos_t){
  const int Np,Ne=5;
  int i;
  FILE *eosfile;
  double T,e,p,cs2,logs;
  eosfile = fopen(filename,"r");
  
  if(eosfile==NULL)
    return 1;
    
  for(i=0;i<Nt;i+=1){
    fscanf(eosfile,"%lf %lf %lf %lf %lf",&T,&cs2,&e,&p,&logs);
    eos_t[i*Ne+0] = logs; eos_t[i*Ne+1] = T;
    eos_t[i*Ne+2] = p;    eos_t[i*Ne+3] = e;
    eos_t[i*Ne+4] = cs2;
  }
  
  fclose(eosfile);
  return 0;
}

int main(){
}
