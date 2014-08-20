#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TCanvas.h"
#include "TEllipse.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TH1D.h"
#include "TLine.h"
#include "TApplication.h"
#else
class TCanvas;
class TEllipse;
class TRandom3;
class TMath;
class TH1D;
class TLine;
class TApplication;
#endif

using namespace std;

class Point {

  protected:
    double fX;
    double fY;

  public:
    Point(double x,double y) : fX(x), fY(y) {  }
    ~Point() { }

    double x() { return fX; }
    double y() { return fY; }

};


class Triangle {

  protected:
    Point* fP1;
    Point* fP2;
    Point* fP3;
    double fAngle;
    double fSize;
    TRandom3* fRand;

  public:
    Triangle() : fAngle(0.), fSize(0.) { fP1 = 0; fP2 = 0; fP3 = 0; fRand = new TRandom3(); fRand->SetSeed(0); }
    ~Triangle() { delete fRand; fRand = 0; }

    void random(double xmin,double xmax,double ymin,double ymax);
    void draw(int color=1);
    bool is_inside(TEllipse* ell);
    double get_size() { return fSize; }
    double get_angle() { return fAngle; }

};

void Triangle::random(double xmin,double xmax,double ymin,double ymax) {

  double x = fRand->Uniform(xmin,xmax);
  double y = fRand->Uniform(ymin,ymax);
  double r = fRand->Uniform(0.3,(xmax-xmin)/2.);
  double a = fRand->Uniform(0.,TMath::TwoPi());

  fSize = r;
  fAngle = a;

  double x1 = x + r*cos(a);
  double y1 = y + r*sin(a);
  if(fP1) delete fP1;
  fP1 = new Point(x1,y1);

  double x2 = x + r*cos(a+1.*TMath::TwoPi()/3.);
  double y2 = y + r*sin(a+1.*TMath::TwoPi()/3.);
  if(fP2) delete fP2;
  fP2 = new Point(x2,y2);

  double x3 = x + r*cos(a+2.*TMath::TwoPi()/3.);
  double y3 = y + r*sin(a+2.*TMath::TwoPi()/3.);
  if(fP3) delete fP3;
  fP3 = new Point(x3,y3);

}

void Triangle::draw(int color) {

  TLine line;
  line.SetLineWidth(2);
  line.SetLineColor(color);
  double x1 = fP1->x();
  double x2 = fP2->x();
  double x3 = fP3->x();
  double y1 = fP1->y();
  double y2 = fP2->y();
  double y3 = fP3->y();
  line.DrawLine(x1,y1,x2,y2);
  line.DrawLine(x1,y1,x3,y3);
  line.DrawLine(x2,y2,x3,y3);

}

bool Triangle::is_inside(TEllipse* ell) {

  double x[3],y[3];
  x[0] = fP1->x();
  x[1] = fP2->x();
  x[2] = fP3->x();
  y[0] = fP1->y();
  y[1] = fP2->y();
  y[2] = fP3->y();

  double x0 = ell->GetX1();
  double y0 = ell->GetY1();
  double rx = ell->GetR1();
  double ry = ell->GetR2();

  int dist = 0;
  for(int i=0;i<3;++i) {
    double contour = pow((x[i]-x0)/rx,2.) + pow((y[i]-y0)/ry,2.);
    if(contour>1.) dist++;
  }

  return dist==0 ? true : false;

}


//*****************************************************************************

void triangles(int ntmax) {

  TCanvas* cc = new TCanvas("cc","Triangles",10,10,600,600);

  TEllipse* ell = new TEllipse(0.5,0.5,0.2,0.4,0,360,0);
  ell->SetFillColor(0);
  ell->SetFillStyle(0);
  ell->SetLineStyle(1);
  ell->SetLineWidth(2);
  ell->Draw();

  cc->Update();

  double xmin = 0.3;
  double xmax = 0.7;
  double ymin = 0.1;
  double ymax = 0.9;

  TH1D* h1 = new TH1D("h1","angular distribution",100,0.,TMath::TwoPi());

  Triangle* tri = new Triangle();
  int nt = 0;
  while(nt<ntmax) {
    tri->random(xmin,xmax,ymin,ymax);
    if(tri->is_inside(ell)) { 
      double angle = tri->get_angle();
      double size = tri->get_size();
      h1->Fill(angle,size*size);
      if(nt<20) tri->draw(2); 
      nt++; 
    }
    //else tri->draw(4);
  }

  TCanvas* c1 = new TCanvas("c1","Correlation",650,10,600,600);
  h1->Draw();

}

///////////////////////////////////////////////////////////////////////////////

#ifndef __CINT__

int main(int argc,char* argv[]) {

  // check arguments
  int ntmax;
  if(argc==1) ntmax = 20;
  if(argc==2) ntmax = atoi(argv[1]);
  if(argc>2) { cout << "Error!" << endl; return -1; }

  //gROOT->SetBatch();
  TApplication *app = new TApplication("app",&argc,argv);

  triangles(ntmax);

  app->Run(true);

  return 0;

}

#endif
