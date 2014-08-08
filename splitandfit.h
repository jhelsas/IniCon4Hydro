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
 
class domain{
  public:
    int D,Nv,type,good;
    vector <double> xv; // eu não deveria precisar de vector aqui
    double S;
  
    int init(int,int,int);
    domain();
    domain(int,int,int);
};

int domain_count(vector <domain> dom);

int domain_check(vector <domain> dom);

int el_print(domain el);

int rollin(int *ind,int *indx,int D,int n);

int rollout(int *ind,int *indx,int D,int n);

int cubic_split(int Nsplit,vector <domain> & dom, int n,int D);

int init_cube(double xl[],double xu[],vector <domain> & dom,int D);

int clean_domain(vector <domain> &dom);

int print_sph(int D,const char *filename,vector <domain> dom);

int domain_split(int D, int Nsplit,double cutoff,vector <domain>& dom, gsl_monte_function F);
