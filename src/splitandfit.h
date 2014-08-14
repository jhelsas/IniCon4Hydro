/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;
 
/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */
 
class domain{
  public:
    int D,Nv,type,good;
    vector <double> xv; 
    //static vector<domain> domains; // -> to be used in a future update
    double S;
  
    int init(int,int,int);
    domain();
    domain(int,int,int);
};

typedef struct wparams{
  domain *mdel;
  void *p;
} wparams;


/*******************************************************
 * Domain Related Functions  
 *
 * Public Functions: to use as 
 *   interface of the library
 * 
 * Internal Functions: to use internally and not
 *   should be called outside the library, 
 *   except for debbuging purposes
 * 
 ********************************************************/

/*
 * Public Functions
 */
 
int domain_count(vector <domain> dom);

int domain_check(vector <domain> dom);

int el_print(domain el);

int print_domain_n(int D,vector <domain> dom,int n);

int print_domain(int D,vector <domain> dom);

int check_inside(double *x,size_t D,domain *mdel);

int init_cube(double xl[],double xu[],vector <domain> & dom,int D);

int init_triangle(int D,double Ntri,double *xv,vector <domain> & dom);

int domain_split(int D,double cutoff,vector <domain>& dom, gsl_monte_function F);

int clean_domain(vector <domain> &dom);

int print_sph(int D,const char *filename,vector <domain> &dom);

/*
 * Internal Functions
 */

int rollin(int *ind,int *indx,int D,int n);

int rollout(int *ind,int *indx,int D,int n);

double aspect_ratio(domain mdel);

int cubic_split(int D,vector <domain> & dom, int n);

int bc_simplex_split(int D,vector <domain> & dom,int n);

int bc_coord(int D,double *r,double *lmb,domain mdel);

int bc_check(int D,double *lmb);

int bc_test(int D,double *r,domain mdel);

int triangle_midpoint_split(int D,vector <domain> & dom,int n);

/**************************************************************
 *  Utility Functions  
 * 
 *  Usefull but not essential for the functioning of the
 *   library. Mainly to facilitatethe init functions.
 * 
 **************************************************************/
 
int unit_hexagon(int Ntri,int D,double *xv);

int unit2_hexagon(int Ntri,int D,double *xv);
