#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include <vector>
#include "../splitandfit.h"

using namespace std;

/*
 * type = 0 : D dimensional cubic domain
 * type = 1 : 2 dimensional triangular domain
 * type = 2 : 3 dimensional tetraedric domain
 */
 
/*
 * Exemplo do uso da Biblioteca de particionamento de domínio.
 * 
 * Os inputs são:
 * 1)Uma função, a ser usada no molde da gsl, 
 *   como as presentes no arquivo "trial-functions.h" ou abaixo tkgauss
 *   Esta função seria posta como membro de uma gsl_monte_function, 
 *   que tem que estar declarada junto com o código.
 * 
 * 2)Uma variavel do tipo wparams, que servirá para compatibilizar a 
 *   função com o domínio de integração. Caso haja algum parâmetro
 *   setável dentro da função, utilize o membro p de wparams para
 *   armazenar os parâmetros, e não o membro params da 
 *   gsl_monte_function, este deve ser reservado para a wparams
 * 
 * 3)Um vetor de domínios dom, o qual será inicializado para ser um 
 *   domínio ou de cubos ou de triangulos, no presente momento
 * 
 * 4)Um valor de cutoff que estabelecará o maior peso (i.e. integral 
 *   da função dada) que um domínio pode ter. Se um dado domínio 
 *   passar deste valor de peso, ele é particionado pelo código.
 * 
 * Todas as funções a serem executadas retornam um int, que é 0 se 
 * a função foi corretamente executada, e 1 caso contrário. Os erros
 * ainda não foram catalogados, mas 0 é sucesso de execução, caso
 * não seja 0, procure no código da função o que representa o erro.
 * 
 * A inicialização por cubos exigem os intervalos aonde está definido
 * O cubo, o que quer dizer:
 * 
 *       xu[1]  _______________________
 *             |                       |
 *             |                       |
 *             |                       |
 *             |                       |
 *             |                       |
 *       xl[1] |_______________________|
 *            xl[0]                   xu[0]
 * 
 * A inicialização por triangulos exige um array de triangulos, com
 * as posições de cada triangulo.
 * 
 * Um dado triangulo é caracterizado por
 * 
 *                 (x[2],x[3])
 *                      /\
 *                     /  \
 *                    /    \
 *                   /      \
 *                  /        \
 *                 /          \
 *                /            \
 *               /______________\
 *          (x[0],x[1])     (x[4],x[5])
 * 
 *  
 * Quando há vários triângulos, numere cada triangulo, e armazene
 * sequencialmente os vertices de cada um deles separadamente, mesmo
 * que haja vértices em comum, isto não é um problema.
 *                  _______________  ______________
 *                 /\              /\              /\
 *                /  \            /  \            /  \
 *               /    \    1     /    \    3     /    \
 *              /      \        /      \        /      \
 *             /   0    \      /   2    \      /    4   \
 *            /          \    /          \    /          \
 *           /            \  /            \  /            \
 *          /______________\/______________\/______________\
 * 
 * Portanto, o elemento do array que dá, em d dimensões 
 * (i.e., há d+1 vertices), a coordenada l do vertice v do triangulo t
 * é:
 * 
 *  x[t*(d*(d+1)) + d*v + l]
 *  
 * Exemplos estão feitos na função unit_hexagon e unit2_hexagon 
 * 
 * Depois disto, as funções domain_split, clean_domain e print_sph
 * devem ser chamadas nesta ordem.
 * 
 * Se clean_domain for chamado antes de domain_split, ela apenas não
 * terá efeito algum, mas se print_sph for chamado antes de clean_domain
 * , ela imprimirá apenas os domínios de input, enquanto se ela for
 * chamada antes de clean_domain mas depois de domain_split, ela irá
 * imprimir partículas SPH vindas de domínios espúrios.
 * 
 * Numa revisão posterior, clean_domain deverá ser imbutida dentro de 
 * domain_split
 */

double cubic_dome(double *x,size_t n,void *par){
  unsigned int i;
  int err;
  double r=0.;
  wparams *lpar=(wparams*)par;
  
  err=check_inside(x,n,lpar->mdel);
  if(err!=0)
    return 0.;
    
  for(i=0;i<=n;i+=1)
    r+=x[i]*x[i];
  r=sqrt(r);
  if(r<=1. && n==1)
    return 0.75*(1.-r*r);
  else if(r<=1. && n==2)
    return (3.*M_1_PI)*(1.-r*r);
  else if(r<=1. && n==3)
    return (1.875*M_1_PI)*(1.-r*r);
  else
    return 0;
  
}

int null_velocity(double *x,size_t dim,void *par,double *u){
  int l;
  for(l=0;l<dim;l+=1)
    u[l]=0.;
  
  return 0;
}

int main(){
  int D=2,Ntri=6,split_type=0;
  int l,err,Npoints,N;
  double cutoff=0.001,xi[D],xf[D],xv[Ntri*(D+1)*D];
  double xl[D],xu[D],dx[D];
  double *xp,*x,*u,*S,s,dist,h=0.1;
  wparams par;
  vector <domain> dom;
  gsl_monte_function F;
  ifstream sphfile;
  ofstream plotfile;
  FILE *sphofile;
  
  F.f = &cubic_dome; 
  F.dim=D; 
  par.p=NULL;
  F.params=&par;
  
  if(split_type==0){
    cout << "init\n";
    for(l=0;l<D;l+=1){xi[l]=-1.;xf[l]=1.;}  
    err=init_cube(xi,xf,dom,D);if(err!=0) return err; 
  }
  else if(split_type==1){
    cout << "unit\n";
    err=unit2_hexagon(Ntri,D,xv);if(err!=0) return err;
    
    cout << "init\n";
    err=init_triangle(D,Ntri,xv,dom);if(err!=0) return err;
  } 
  
  cout << "split\n";
  err=domain_split(D,cutoff,dom,F); if(err!=0){ cout << "out: " << err << endl;return err;}
  
  cout << "clean\n";
  err=clean_domain(dom);if(err!=0) return err;
  
  cout << "print\n";  
  err=print_moving_sph(D,"results/cubic_dome.dat",dom,null_velocity,&par); if(err!=0) return err;
    
  cout << "reading\n";
  err=sph_read("results/cubic_dome.dat",&D,&N,&x,&u,&S);if(err!=0) return err;

  for(l=0;l<D;l+=1){xl[l]=-4.0;dx[l]=0.15;xu[l]=4.0+1.01*dx[l];}
  
  err=create_grid(D,&xp,xl,xu,dx,&Npoints);if(err!=0) return err;         
  
  cout << "ploting\n";
  err=sph_dens(D,N,Npoints,xp,x,S,h,xl,xu,cubic_dome,"results/cubic_dome_plot.dat",NULL);if(err!=0){ cout << err << endl; return err;}
  
  delete x;delete u; delete S; 
  
  return 0;
}
