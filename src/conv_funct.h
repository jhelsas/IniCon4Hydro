typedef struct conv_wrap{
  double (*f2s)(double,void*p);
  void *funpar, *f2spar;
} conv_wrap;

double e2s_pion(double epsilon,void *p);

double e2s_qg(double epsilon,void *p);

double e2s_qgphr(double epsilon, void *p);

double e2s_table(double epsilon,void *p);

double T2s_table(double T,void *p);
