#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

typedef struct conv_wrap{
  double (*f2s)(double,void*p);
  void *funpar, *f2spar;
} conv_wrap;

typedef struct eosp{
  double T,s,p,e,h,hsh,cs2;
} eosp;

extern gsl_interp_accel *acc; // global variables
extern gsl_spline *splT, *splp ,*sple ,*splcs2 ;
extern const gsl_root_fsolver_type *rsT;
extern gsl_root_fsolver *solver;
extern gsl_interp_accel *global_e2s_zoltan_acc; // global variables
extern gsl_spline *global_e2s_zoltan_spls;

double e2s_pion(double epsilon,void *p);

double e2s_qg(double epsilon,void *p);

double e2s_qgphr(double epsilon, void *p);

double e2s_zoltan(double e, void *p);

double e2s_table(double epsilon,void *p);

double T2s_table(double T,void *p);
