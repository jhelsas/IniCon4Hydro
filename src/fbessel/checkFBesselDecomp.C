//#define DEBUG
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TObjectTable.h"
#include "FourierBesselDecomp2D.h"
#else
class TROOT;
class TApplication;
class TSystem;
class TStyle;
class TFile;
class TChain;
class TCanvas;
class TString;
class TH1;
class TH2;
class TH3;
class TGraph;
class TProfile2D;
class TTree;
class TColor;
class TGaxis;
class TObjectTable;
class FourierBesselDecomp2D;
#endif

using namespace::std;

void checkFBesselDecomp(int mMin=-8,int mMax=8,     // range of Jm(x) to be included in the decomposition
                        int nMin=1, int nMax=8) {   // range of Jm(x) zeros to be computed


#ifdef __CINT__
  // compile/load needed libraries
  cout << "loading libraries..." << endl;
  gSystem->Load("libMathMore");
  gROOT->LoadMacro("FourierBesselDecomp2D.cxx+");
  cout << "libraries loaded!" << endl;
#endif

  gStyle->SetPalette(1);
  gStyle->SetNumberContours(40);
  //  gStyle->SetCanvasPreferGL(true);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetStatColor(0);
  gStyle->SetStatFont(132);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(132,"t");
  gStyle->SetNumberContours(100);

  gStyle->SetPalette(1, 0);
  const int NCont = 50; // should suffice, can be set to 255
  const int NRGBs = 5;
  double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  double red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  double green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  double blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  //
  // USE THIS FOR GRAYSCALE PLOTS
  //const int NCont = 25; // should suffice, can be set to 255
  //const int NRGBs = 7;
  //double stops[NRGBs] = { 0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00 };
  //double red[NRGBs] = { 1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
  //double green[NRGBs] = { 1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
  //double blue[NRGBs] = { 1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
  //
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  TGaxis::SetMaxDigits(3);

  // define canvases to draw the results
  TCanvas* c1 = new TCanvas("c1","",70,15,350,350);
  c1->SetTicks();
  c1->SetLeftMargin(0.20);
  c1->SetBottomMargin(0.20);
  c1->SetRightMargin(0.16);
  c1->SetTopMargin(0.10);
  c1->SetGridx();
  c1->SetGridy();

  TCanvas* c2 = new TCanvas("c2","",430,15,350,350);
  c2->SetTicks();
  c2->SetLeftMargin(0.20);
  c2->SetBottomMargin(0.20);
  c2->SetRightMargin(0.16);
  c2->SetTopMargin(0.10);
  c2->SetGridx();
  c2->SetGridy();

  TCanvas* c3 = new TCanvas("c3","",790,15,350,350);
  c3->SetTicks();
  c3->SetLeftMargin(0.20);
  c3->SetBottomMargin(0.20);
  c3->SetRightMargin(0.20);
  c3->SetTopMargin(0.10);

  TCanvas* c4 = new TCanvas("c4","",70,400,350,350);
  c4->SetTicks();
  c4->SetLeftMargin(0.20);
  c4->SetBottomMargin(0.20);
  c4->SetRightMargin(0.05);
  c4->SetTopMargin(0.05);

  TCanvas* c5 = new TCanvas("c5","",430,400,350,350);
  c5->SetTicks();
  c5->SetLeftMargin(0.20);
  c5->SetBottomMargin(0.20);
  c5->SetRightMargin(0.05);
  c5->SetTopMargin(0.05);

  TCanvas* c6 = new TCanvas("c6","",790,400,350,350);
  c6->SetTicks();
  c6->SetLeftMargin(0.20);
  c6->SetBottomMargin(0.20);
  c6->SetRightMargin(0.05);
  c6->SetTopMargin(0.05);

  // define histogram for the tranverse profile of the initial condition
  TH2D* h2d = new TH2D("h2d","Distribuic#tilde{a}o Original",25,-8.,8.,25,-8.,8.);
  //TH2D* h2d = new TH2D("h2d","Distribuic#tilde{a}o Original",100,-8.,8.,100,-8.,8.);
  h2d->SetStats(false);
  h2d->GetXaxis()->SetTitle("x [fm]");
  h2d->GetXaxis()->SetTitleFont(132);
  h2d->GetXaxis()->SetTitleSize(0.06);
  h2d->GetXaxis()->SetTitleOffset(1.25);
  h2d->GetXaxis()->CenterTitle(true);
  h2d->GetXaxis()->SetLabelFont(132);
  h2d->GetXaxis()->SetLabelSize(0.05);
  h2d->GetXaxis()->SetLabelOffset(0.015);
  h2d->GetXaxis()->SetNdivisions(310);
  h2d->GetYaxis()->SetTitle("y [fm]");
  h2d->GetYaxis()->SetTitleFont(132);
  h2d->GetYaxis()->SetTitleSize(0.06);
  h2d->GetYaxis()->SetTitleOffset(1.35);
  h2d->GetYaxis()->CenterTitle(true);
  h2d->GetYaxis()->SetLabelFont(132);
  h2d->GetYaxis()->SetLabelSize(0.05);
  h2d->GetYaxis()->SetLabelOffset(0.015);
  h2d->GetYaxis()->SetNdivisions(310);
  h2d->GetZaxis()->SetTitleFont(132);
  h2d->GetZaxis()->SetTitleSize(0.06);
  h2d->GetZaxis()->SetTitleOffset(1.35);
  h2d->GetZaxis()->CenterTitle(true);
  h2d->GetZaxis()->SetLabelFont(132);
  h2d->GetZaxis()->SetLabelSize(0.05);
  h2d->GetZaxis()->SetLabelOffset(0.015);
  h2d->GetZaxis()->SetNdivisions(310);

  // define histograms for the Fourier-Bessel norms
  TH1D* hL2 = new TH1D("L2","",1000,0.,1000.);
  TH1D* hH1 = new TH1D("H1","",1000,0.,1000.);
  TH1D* hM1 = new TH1D("M1","",1000,0.,1000.);

  // define graphs for the Fourier-Bessel norms
  TGraph* grL2ang = new TGraph();
  TGraph* grL2rad = new TGraph();
  TGraph* grH1ang = new TGraph();
  TGraph* grH1rad = new TGraph();
  TGraph* grM1ang = new TGraph();
  TGraph* grM1rad = new TGraph();

  // define object to compute 2D Fourier-Bessel decomposition
  FourierBesselDecomp2D* fourier2D = new FourierBesselDecomp2D();
  fourier2D->SetRadius(20.);
  fourier2D->SetAngularCoefsRange(mMin,mMax);
  fourier2D->SetRadialCoefsRange(nMin,nMax);

  // define histogram to store runtime
  TH1D* hRunTime = new TH1D("hRunTime","",10000,0.,10.);

  int nbx = h2d->GetNbinsX();
  int nby = h2d->GetNbinsY();

  //  TF2* fang = new TF2("fang","10.*exp(-0.05*(x*x+y*y))*sin(5.*[0]*atan2(y,x))",-8.,8.,-8.,8.);
  TF2* fang = new TF2("fang","10.*sin(5.*[0]*atan2(y,x))",-8.,8.,-8.,8.);
  fang->SetNpx(nbx);
  fang->SetNpy(nby);
  //  TF2* frad = new TF2("frad","10.*exp(-0.05*(x*x+y*y))*sin([0]*sqrt(x*x+y*y))",-8.,8.,-8.,8.);
  TF2* frad = new TF2("frad","10.*sin([0]*sqrt(x*x+y*y))",-8.,8.,-8.,8.);
  frad->SetNpx(nbx);
  frad->SetNpy(nby);

  // open file to save the plots
  TFile* fout = new TFile("fbesselCheck.root","RECREATE");

  for(int ipar=1;ipar<=10;++ipar) {
    double par = 0.2*double(ipar);

    // First do angular test
    fang->SetParameter(0,par);
    h2d->Reset();
    for(int ibx=0;ibx<nbx;++ibx) {
      for(int iby=0;iby<nby;++iby) {
	double x = h2d->GetXaxis()->GetBinCenter(ibx+1);
	double y = h2d->GetYaxis()->GetBinCenter(iby+1);
	double val = fang->Eval(x,y);
	h2d->Fill(x,y,val);
      }
    }

    // draw original initial condition transverse profile
    c1->cd();
    cout << "  Ploting the resulting profile..." << endl;
    h2d->Draw("colz");
    h2d->Write(Form("histOriginal_ang%d%d",(ipar)/10,(ipar)%10));

    // draw resulting Fourier-Bessel decomposition of the original distribution
    c2->cd();
    cout << "  Now computing the 2D decomposition... (it may take some time)";
    fourier2D->DoDecomposition(h2d);
    cout << "\r  Now computing the 2D decomposition... done!                     " << endl;
    cout << "  Ploting the results..." << endl;
    //    fourier2D->Draw("colz");
    TF2* f2d = (TF2*)fourier2D->GetFunction();
    // TH2D* hframe = (TH2D*)f2d->GetHistogram();
    // hframe->SetStats(false);
    f2d->GetHistogram()->SetTitle("2D Fourier-Bessel Decomposition");
    f2d->GetHistogram()->GetXaxis()->SetTitle("x [fm]");
    f2d->GetHistogram()->GetXaxis()->SetTitleFont(132);
    f2d->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    f2d->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    f2d->GetHistogram()->GetXaxis()->CenterTitle(true);
    f2d->GetHistogram()->GetXaxis()->SetLabelFont(132);
    f2d->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    f2d->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    f2d->GetHistogram()->GetXaxis()->SetNdivisions(310);
    f2d->GetHistogram()->GetYaxis()->SetTitle("y [fm]");
    f2d->GetHistogram()->GetYaxis()->SetTitleFont(132);
    f2d->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    f2d->GetHistogram()->GetYaxis()->SetTitleOffset(1.35);
    f2d->GetHistogram()->GetYaxis()->CenterTitle(true);
    f2d->GetHistogram()->GetYaxis()->SetLabelFont(132);
    f2d->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    f2d->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    f2d->GetHistogram()->GetYaxis()->SetNdivisions(310);
    f2d->GetHistogram()->GetZaxis()->SetTitleFont(132);
    f2d->GetHistogram()->GetZaxis()->SetTitleSize(0.06);
    f2d->GetHistogram()->GetZaxis()->SetTitleOffset(1.35);
    f2d->GetHistogram()->GetZaxis()->CenterTitle(true);
    f2d->GetHistogram()->GetZaxis()->SetLabelFont(132);
    f2d->GetHistogram()->GetZaxis()->SetLabelSize(0.05);
    f2d->GetHistogram()->GetZaxis()->SetLabelOffset(0.015);
    f2d->GetHistogram()->GetZaxis()->SetNdivisions(310);
    f2d->Draw("colz");
    f2d->Write(Form("funcFBdecomp_ang%d%d",(ipar)/10,(ipar)%10));

    c1->Modified();
    c1->Update();
    c1->Write(Form("DistOrig_ang%d%d",(ipar)/10,(ipar)%10));
    c1->SaveAs(Form("DistOrig_ang%d%d.eps",(ipar)/10,(ipar)%10));
    c2->Modified();
    c2->Update();
    c2->Write(Form("DistFBdecomp_ang%d%d",(ipar)/10,(ipar)%10));
    c2->SaveAs(Form("DistFBdecomp_ang%d%d.eps",(ipar)/10,(ipar)%10));

    // draw Bessel Jm(x) zeros
    c3->cd();
    fourier2D->DrawZerosHisto("colz");
    c3->Modified();
    c3->Update();

    // draw real part of the coefficients Amn
    // c4->cd();
    // fourier2D->DrawReCoefsHisto("colz");
    // c4->Modified();
    // c4->Update();

    // draw imaginary part of the coefficients Amn
    // c5->cd();
    // fourier2D->DrawImCoefsHisto("colz");
    // c5->Modified();
    // c5->Update();

    // draw the norm of the coefficients Amn
    // c6->cd();
    // fourier2D->DrawCoefsNormHisto("colz");
    // c6->Modified();
    // c6->Update();

    // fill histograms of the norms
    double l2 = fourier2D->GetL2();
    hL2->Fill(l2);
    grL2ang->SetPoint(ipar-1,par,l2);
  
    double h1 = fourier2D->GetH1();
    hH1->Fill(h1);
    grH1ang->SetPoint(ipar-1,par,h1);
    
    double m1 = fourier2D->GetM1();
    hM1->Fill(m1);
    grM1ang->SetPoint(ipar-1,par,m1);
  
    // draw H1
    c4->cd();
    //c4->Clear();
    //hL2->Draw();
    grL2ang->SetMarkerStyle(24);
    grL2ang->GetHistogram()->GetXaxis()->SetTitle("#omega_{ang} [a.u.]");
    grL2ang->GetHistogram()->GetXaxis()->SetTitleFont(132);
    grL2ang->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    grL2ang->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    grL2ang->GetHistogram()->GetXaxis()->CenterTitle(true);
    grL2ang->GetHistogram()->GetXaxis()->SetLabelFont(132);
    grL2ang->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    grL2ang->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    grL2ang->GetHistogram()->GetXaxis()->SetNdivisions(310);
    grL2ang->GetHistogram()->GetYaxis()->SetTitle("L2");
    grL2ang->GetHistogram()->GetYaxis()->SetTitleFont(132);
    grL2ang->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    grL2ang->GetHistogram()->GetYaxis()->SetTitleOffset(1.60);
    grL2ang->GetHistogram()->GetYaxis()->CenterTitle(true);
    grL2ang->GetHistogram()->GetYaxis()->SetLabelFont(132);
    grL2ang->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    grL2ang->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    grL2ang->GetHistogram()->GetYaxis()->SetNdivisions(310);
    if(ipar==1) grL2ang->Draw("alp");
    else grL2ang->Draw("lp");
    c4->Modified();
    c4->Update();

    // draw M1
    c5->cd();
    //c5->Clear();
    //hH1->Draw();
    grH1ang->SetMarkerStyle(25);
    grH1ang->SetMarkerColor(2);
    grH1ang->SetLineColor(2);
    grH1ang->GetHistogram()->GetXaxis()->SetTitle("#omega_{ang} [a.u.]");
    grH1ang->GetHistogram()->GetXaxis()->SetTitleFont(132);
    grH1ang->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    grH1ang->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    grH1ang->GetHistogram()->GetXaxis()->CenterTitle(true);
    grH1ang->GetHistogram()->GetXaxis()->SetLabelFont(132);
    grH1ang->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    grH1ang->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    grH1ang->GetHistogram()->GetXaxis()->SetNdivisions(310);
    grH1ang->GetHistogram()->GetYaxis()->SetTitle("H1");
    grH1ang->GetHistogram()->GetYaxis()->SetTitleFont(132);
    grH1ang->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    grH1ang->GetHistogram()->GetYaxis()->SetTitleOffset(1.60);
    grH1ang->GetHistogram()->GetYaxis()->CenterTitle(true);
    grH1ang->GetHistogram()->GetYaxis()->SetLabelFont(132);
    grH1ang->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    grH1ang->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    grH1ang->GetHistogram()->GetYaxis()->SetNdivisions(310);
    if(ipar==1) grH1ang->Draw("alp");
    else grH1ang->Draw("lp");
    c5->Modified();
    c5->Update();
    
    // draw M1
    c6->cd();
    //c6->Clear();
    //hM1->Draw();
    grM1ang->SetMarkerStyle(26);
    grM1ang->SetMarkerColor(4);
    grM1ang->SetLineColor(4);
    grM1ang->GetHistogram()->GetXaxis()->SetTitle("#omega_{ang} [a.u.]");
    grM1ang->GetHistogram()->GetXaxis()->SetTitleFont(132);
    grM1ang->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    grM1ang->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    grM1ang->GetHistogram()->GetXaxis()->CenterTitle(true);
    grM1ang->GetHistogram()->GetXaxis()->SetLabelFont(132);
    grM1ang->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    grM1ang->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    grM1ang->GetHistogram()->GetXaxis()->SetNdivisions(310);
    grM1ang->GetHistogram()->GetYaxis()->SetTitle("M1");
    grM1ang->GetHistogram()->GetYaxis()->SetTitleFont(132);
    grM1ang->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    grM1ang->GetHistogram()->GetYaxis()->SetTitleOffset(1.60);
    grM1ang->GetHistogram()->GetYaxis()->CenterTitle(true);
    grM1ang->GetHistogram()->GetYaxis()->SetLabelFont(132);
    grM1ang->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    grM1ang->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    grM1ang->GetHistogram()->GetYaxis()->SetNdivisions(310);
    if(ipar==1) grM1ang->Draw("alp");
    else grM1ang->Draw("lp");
    c6->Modified();
    c6->Update();

  }

  c4->Write("normL2_angFluct");
  c5->Write("normH1_angFluct");
  c6->Write("normM1_angFluct");

  grL2ang->Write("grNormL2_angFluct");
  grH1ang->Write("grNormH1_angFluct");
  grM1ang->Write("grNormM1_angFluct");

  c4->SaveAs(Form("normL2_angFluct.eps"));
  c5->SaveAs(Form("normH1_angFluct.eps"));
  c6->SaveAs(Form("normM1_angFluct.eps"));

  c1->Clear();
  c2->Clear();
  c3->Clear();
  c4->Clear();
  c5->Clear();
  c6->Clear();

  for(int ipar=1;ipar<=10;++ipar) {
    double par = 0.2*double(ipar);

    // Now do the radial part
    frad->SetParameter(0,par);
    h2d->Reset();
    for(int ibx=0;ibx<nbx;++ibx) {
      for(int iby=0;iby<nby;++iby) {
	double x = h2d->GetXaxis()->GetBinCenter(ibx+1);
	double y = h2d->GetYaxis()->GetBinCenter(iby+1);
	double val = frad->Eval(x,y);
	h2d->Fill(x,y,val);
      }
    }

    // draw original initial condition transverse profile
    c1->cd();
    cout << "  Ploting the resulting profile..." << endl;
    h2d->Draw("colz");
    h2d->Write(Form("histOriginal_rad%d%d",(ipar)/10,(ipar)%10));

    // draw resulting Fourier-Bessel decomposition of the original distribution
    c2->cd();
    cout << "  Now computing the 2D decomposition... (it may take some time)";
    fourier2D->DoDecomposition(h2d);
    cout << "\r  Now computing the 2D decomposition... done!                     " << endl;
    cout << "  Ploting the results..." << endl;
    //    fourier2D->Draw("colz");
    TF2* f2d = (TF2*)fourier2D->GetFunction();
    // TH2D* hframe = (TH2D*)f2d->GetHistogram();
    // hframe->SetStats(false);
    f2d->GetHistogram()->SetTitle("2D Fourier-Bessel Decomposition");
    f2d->GetHistogram()->GetXaxis()->SetTitle("x [fm]");
    f2d->GetHistogram()->GetXaxis()->SetTitleFont(132);
    f2d->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    f2d->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    f2d->GetHistogram()->GetXaxis()->CenterTitle(true);
    f2d->GetHistogram()->GetXaxis()->SetLabelFont(132);
    f2d->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    f2d->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    f2d->GetHistogram()->GetXaxis()->SetNdivisions(310);
    f2d->GetHistogram()->GetYaxis()->SetTitle("y [fm]");
    f2d->GetHistogram()->GetYaxis()->SetTitleFont(132);
    f2d->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    f2d->GetHistogram()->GetYaxis()->SetTitleOffset(1.35);
    f2d->GetHistogram()->GetYaxis()->CenterTitle(true);
    f2d->GetHistogram()->GetYaxis()->SetLabelFont(132);
    f2d->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    f2d->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    f2d->GetHistogram()->GetYaxis()->SetNdivisions(310);
    f2d->GetHistogram()->GetZaxis()->SetTitleFont(132);
    f2d->GetHistogram()->GetZaxis()->SetTitleSize(0.06);
    f2d->GetHistogram()->GetZaxis()->SetTitleOffset(1.35);
    f2d->GetHistogram()->GetZaxis()->CenterTitle(true);
    f2d->GetHistogram()->GetZaxis()->SetLabelFont(132);
    f2d->GetHistogram()->GetZaxis()->SetLabelSize(0.05);
    f2d->GetHistogram()->GetZaxis()->SetLabelOffset(0.015);
    f2d->GetHistogram()->GetZaxis()->SetNdivisions(310);
    f2d->Draw("colz");
    f2d->Write(Form("funcFBdecomp_rad%d%d",(ipar)/10,(ipar)%10));
    c1->Modified();
    c1->Update();
    c1->Write(Form("DistOrig_rad%d%d",(ipar)/10,(ipar)%10));
    c1->SaveAs(Form("DistOrig_rad%d%d.eps",(ipar)/10,(ipar)%10));
    c2->Modified();
    c2->Update();
    c2->Write(Form("DistFBdecomp_rad%d%d",(ipar)/10,(ipar)%10));
    c2->SaveAs(Form("DistFBdecomp_rad%d%d.eps",(ipar)/10,(ipar)%10));
    
    // draw Bessel Jm(x) zeros
    c3->cd();
    fourier2D->DrawZerosHisto("colz");
    c3->Modified();
    c3->Update();

    // draw real part of the coefficients Amn
    // c4->cd();
    // fourier2D->DrawReCoefsHisto("colz");
    // c4->Modified();
    // c4->Update();

    // draw imaginary part of the coefficients Amn
    // c5->cd();
    // fourier2D->DrawImCoefsHisto("colz");
    // c5->Modified();
    // c5->Update();

    // draw the norm of the coefficients Amn
    // c6->cd();
    // fourier2D->DrawCoefsNormHisto("colz");
    // c6->Modified();
    // c6->Update();

    // fill histograms of the norms
    double l2 = fourier2D->GetL2();
    hL2->Fill(l2);
    grL2rad->SetPoint(ipar-1,par,l2);
  
    double h1 = fourier2D->GetH1();
    hH1->Fill(h1);
    grH1rad->SetPoint(ipar-1,par,h1);
    
    double m1 = fourier2D->GetM1();
    hM1->Fill(m1);
    grM1rad->SetPoint(ipar-1,par,m1);
  
    // draw H1
    c4->cd();
    //c4->Clear();
    //hL2->Draw();
    grL2rad->SetMarkerStyle(24);
    grL2rad->GetHistogram()->GetXaxis()->SetTitle("#omega_{rad} [a.u.]");
    grL2rad->GetHistogram()->GetXaxis()->SetTitleFont(132);
    grL2rad->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    grL2rad->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    grL2rad->GetHistogram()->GetXaxis()->CenterTitle(true);
    grL2rad->GetHistogram()->GetXaxis()->SetLabelFont(132);
    grL2rad->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    grL2rad->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    grL2rad->GetHistogram()->GetXaxis()->SetNdivisions(310);
    grL2rad->GetHistogram()->GetYaxis()->SetTitle("L2");
    grL2rad->GetHistogram()->GetYaxis()->SetTitleFont(132);
    grL2rad->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    grL2rad->GetHistogram()->GetYaxis()->SetTitleOffset(1.60);
    grL2rad->GetHistogram()->GetYaxis()->CenterTitle(true);
    grL2rad->GetHistogram()->GetYaxis()->SetLabelFont(132);
    grL2rad->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    grL2rad->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    grL2rad->GetHistogram()->GetYaxis()->SetNdivisions(310);
    if(ipar==1) grL2rad->Draw("alp");
    else grL2rad->Draw("lp");
    c4->Modified();
    c4->Update();

    // draw M1
    c5->cd();
    //c5->Clear();
    //hH1->Draw();
    grH1rad->SetMarkerStyle(25);
    grH1rad->SetMarkerColor(2);
    grH1rad->SetLineColor(2);
    grH1rad->GetHistogram()->GetXaxis()->SetTitle("#omega_{rad} [a.u.]");
    grH1rad->GetHistogram()->GetXaxis()->SetTitleFont(132);
    grH1rad->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    grH1rad->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    grH1rad->GetHistogram()->GetXaxis()->CenterTitle(true);
    grH1rad->GetHistogram()->GetXaxis()->SetLabelFont(132);
    grH1rad->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    grH1rad->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    grH1rad->GetHistogram()->GetXaxis()->SetNdivisions(310);
    grH1rad->GetHistogram()->GetYaxis()->SetTitle("H1");
    grH1rad->GetHistogram()->GetYaxis()->SetTitleFont(132);
    grH1rad->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    grH1rad->GetHistogram()->GetYaxis()->SetTitleOffset(1.60);
    grH1rad->GetHistogram()->GetYaxis()->CenterTitle(true);
    grH1rad->GetHistogram()->GetYaxis()->SetLabelFont(132);
    grH1rad->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    grH1rad->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    grH1rad->GetHistogram()->GetYaxis()->SetNdivisions(310);
    if(ipar==1) grH1rad->Draw("alp");
    else grH1rad->Draw("lp");
    c5->Modified();
    c5->Update();
    
    // draw M1
    c6->cd();
    //c6->Clear();
    //hM1->Draw();
    grM1rad->SetMarkerStyle(26);
    grM1rad->SetMarkerColor(4);
    grM1rad->SetLineColor(4);
    grM1rad->GetHistogram()->GetXaxis()->SetTitle("#omega_{rad} [a.u.]");
    grM1rad->GetHistogram()->GetXaxis()->SetTitleFont(132);
    grM1rad->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
    grM1rad->GetHistogram()->GetXaxis()->SetTitleOffset(1.25);
    grM1rad->GetHistogram()->GetXaxis()->CenterTitle(true);
    grM1rad->GetHistogram()->GetXaxis()->SetLabelFont(132);
    grM1rad->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
    grM1rad->GetHistogram()->GetXaxis()->SetLabelOffset(0.015);
    grM1rad->GetHistogram()->GetXaxis()->SetNdivisions(310);
    grM1rad->GetHistogram()->GetYaxis()->SetTitle("M1");
    grM1rad->GetHistogram()->GetYaxis()->SetTitleFont(132);
    grM1rad->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
    grM1rad->GetHistogram()->GetYaxis()->SetTitleOffset(1.60);
    grM1rad->GetHistogram()->GetYaxis()->CenterTitle(true);
    grM1rad->GetHistogram()->GetYaxis()->SetLabelFont(132);
    grM1rad->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
    grM1rad->GetHistogram()->GetYaxis()->SetLabelOffset(0.015);
    grM1rad->GetHistogram()->GetYaxis()->SetNdivisions(310);
    if(ipar==1) grM1rad->Draw("alp");
    else grM1rad->Draw("lp");
    c6->Modified();
    c6->Update();

  }

  c4->Write("normL2_radFluct");
  c5->Write("normH1_radFluct");
  c6->Write("normM1_radFluct");

  grL2rad->Write("grNormL2_radFluct");
  grH1rad->Write("grNormH1_radFluct");
  grM1rad->Write("grNormM1_radFluct");

  c4->SaveAs(Form("normL2_radFluct.eps"));
  c5->SaveAs(Form("normH1_radFluct.eps"));
  c6->SaveAs(Form("normM1_radFluct.eps"));

  //delete fourier2D; fourier2D = 0;

}


///////////////////////////////////////////////////////////////////////////////

void printUsage() {
  printf("\nUsage:\n\n");
  printf("      doICoFourierBesselDecomp [options] <out.root file> <sphEventTree_1.root> <sphEventTree_2.root> ... <sphEventTree_N.root> \n");
  printf("      doICoFourierBesselDecomp [options] <out.root file> <inputFiles.list> \n");
  printf("      doICoFourierBesselDecomp [options] <inputFiles.list> \n\n");
  printf("      options:\n");
  printf("      --draw    draw plots during the loop over the events\n");
  printf("      --png     save png files of the plots\n\n");
}


int main(int argc, char *argv[]) {

  TApplication *app = new TApplication("app",&argc,argv);

  checkFBesselDecomp();

  app->Run(true);

  return 0;

}

