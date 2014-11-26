
int ini_write(int D,const char *filename, char *hash, char *folder,
              vector <domain> &dom,
              int (*velocity)(double *,size_t,void *,double *),
              void *par, config *cfg);
