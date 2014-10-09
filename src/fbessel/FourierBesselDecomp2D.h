// Author: Rafael D. de Souza 31/05/2012 

#ifndef FOURIERBESSELDECOMP2D_H
#define FOURIERBESSELDECOMP2D_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// FourierBesselDecomp2D.h - Class definition of FourierBesselDecomp2D  //
//                                                                      //
// Description:                                                         //
// ____________________________________________________________________ //
//                                                                      //
// Last modified: 31/01/2013 - Rafael D. de Souza                       //
// - added setters for m_min m_max n_min n_max r0                       //
//                                                                      //
// Last modified: 31/05/2012 - Rafael D. de Souza                       //
// - first version release                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TString.h"
#include "TMath.h"
#include "Math/SpecFunc.h"
#include "BesselJmZeros.h"

using namespace::std;

class FourierBesselDecomp2D {

 protected:
  
  int m_min;
  int m_max;
  int n_min;
  int n_max;
  double r0;
  
  double BesselZero(int m,int n);
  double Re_phi_mn(int m,int n,double x,double y);
  double Im_phi_mn(int m,int n,double x,double y);
  double ComputeReAmn(int m,int n,TH2D* h2d);
  double ComputeImAmn(int m,int n,TH2D* h2d);
  double func(double* x,double* Amn);
  TGraph2D* ConvertToGraph2D(TH2D* h2d);
  
  TF2* f2d;
  TH2D* hZeros;
  TH2D* hReCoefs;
  TH2D* hImCoefs;
  TH2D* hCoefsNorm;
   
  
 public:
  
  FourierBesselDecomp2D();
  virtual ~FourierBesselDecomp2D();
  
  void DoDecomposition(TH2D* h2d);
  void Draw(TString opt="");
  void DrawZerosHisto(TString opt="");
  void DrawReCoefsHisto(TString opt="");
  void DrawImCoefsHisto(TString opt="");
  void DrawCoefsNormHisto(TString opt="");
  void Plot_phi_mn(int m,int n,double xmin,double xmax,double ymin,double ymax,TString comp="Re");
  void PrintCoefs();

  void SetRadius(double radius);
  void SetAngularCoefsRange(double mMin,double mMax);
  void SetRadialCoefsRange(double nMin,double nMax);
  
  double GetL2();
  double GetH1();
  double GetM1();

  TF2* GetFunction();

  ClassDef(FourierBesselDecomp2D,0)  // FourierBesselDecomp2D class

};

inline void FourierBesselDecomp2D::PrintCoefs() { f2d->Print(); }
inline void FourierBesselDecomp2D::SetRadius(double radius) { r0 = radius; }
inline void FourierBesselDecomp2D::SetAngularCoefsRange(double mMin,double mMax) { m_min = mMin; m_max = mMax; }
inline void FourierBesselDecomp2D::SetRadialCoefsRange(double nMin,double nMax) { n_min = nMin; n_max = nMax; }
inline TF2* FourierBesselDecomp2D::GetFunction() { return (TF2*)f2d; }
  
#endif
