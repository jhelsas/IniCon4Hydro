// Author: Rafael D. de Souza   31/05/2012

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// FourierBesselDecomp2D.cxx - Class implementation of FourierBesselDecomp2D //
//                                                                           //
// Author: Rafael D. de Souza 31/05/2012 (rderradi@ifi.unicamp.br)           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "FourierBesselDecomp2D.h"

ClassImp(FourierBesselDecomp2D)


//______________________________________________________________________________
FourierBesselDecomp2D::FourierBesselDecomp2D() : m_min(-8),m_max(8),n_min(1),n_max(8),r0(12.)
{
  // default constructor
  
  f2d = 0;
  hZeros = 0;
  hReCoefs = 0;
  hImCoefs = 0;
  hCoefsNorm = 0;
  
}


//______________________________________________________________________________
FourierBesselDecomp2D::~FourierBesselDecomp2D()
{
  // default destructor
  
  if(f2d)        { delete f2d;        f2d = 0; }
  if(hZeros)     { delete hZeros;     hZeros = 0; }
  if(hReCoefs)   { delete hReCoefs;   hReCoefs = 0; }
  if(hImCoefs)   { delete hImCoefs;   hImCoefs = 0; }
  if(hCoefsNorm) { delete hCoefsNorm; hCoefsNorm = 0; }
  
}


//______________________________________________________________________________
double FourierBesselDecomp2D::BesselZero(int m,int n)
{
  // return the nth zero of Jm(x)
  
  double lamb_mn = lambda[TMath::Abs(m)][n]; // from the table BesselJmZeros.h
  
  return lamb_mn;
  
}


//______________________________________________________________________________
double FourierBesselDecomp2D::Re_phi_mn(int m,int n,double x,double y)
{
  // real part of phi_mn(x,y)
  
  double r = TMath::Sqrt(pow(x,2.)+pow(y,2.)); r /= r0; // using unit radius
  double theta = TMath::ATan2(y,x);
  double lamb_mn = BesselZero(m,n);
  double _m = (double)m;
  double _n = (double)n;
  double phi_mn = (pow(_m/TMath::Abs(_m),TMath::Abs(_m)))*(ROOT::Math::cyl_bessel_j(TMath::Abs(_m),r*lamb_mn))/(ROOT::Math::cyl_bessel_j(TMath::Abs(_m)+1.,lamb_mn))*(TMath::Cos(_m*theta));
  hZeros->Fill(_m,_n,lamb_mn);

  return phi_mn;
}


//______________________________________________________________________________
double FourierBesselDecomp2D::Im_phi_mn(int m,int n,double x,double y)
{
  // imaginary part of phi_mn(x,y)
  
  double r = TMath::Sqrt(pow(x,2.)+pow(y,2.)); r /= r0; // using unit radius
  double theta = TMath::ATan2(y,x);
  double lamb_mn = BesselZero(m,n);
  double _m = (double)m;
  double _n = (double)n;
  double phi_mn = (pow(_m/TMath::Abs(_m),TMath::Abs(_m)))*(ROOT::Math::cyl_bessel_j(TMath::Abs(_m),r*lamb_mn))/(ROOT::Math::cyl_bessel_j(TMath::Abs(_m)+1.,lamb_mn))*(TMath::Sin(_m*theta));
  hZeros->Fill(_m,_n,lamb_mn);
  
  return phi_mn;
}


//______________________________________________________________________________
double FourierBesselDecomp2D::ComputeReAmn(int m,int n,TH2D* h2d)
{
  // compute Re{Amn} coeficients

  double Amn = 0.;
  int nix = h2d->GetNbinsX();
  int niy = h2d->GetNbinsY();
  double xcm = 0.;//h2d->GetMean(1);
  double ycm = 0.;//h2d->GetMean(2);
  double dx = h2d->GetXaxis()->GetBinWidth(1);
  double dy = h2d->GetYaxis()->GetBinWidth(1);
  for(int ix=0;ix<nix;++ix) {
    for(int iy=0;iy<niy;++iy) {
      double x = h2d->GetXaxis()->GetBinCenter(ix+1) - xcm;
      double y = h2d->GetYaxis()->GetBinCenter(iy+1) - ycm;
      double fxy = h2d->GetBinContent(ix+1,iy+1);
      Amn += (1./(TMath::Pi()*r0*r0))*(fxy*Re_phi_mn(m,n,x,y)*dx*dy);
    }
  }
  
  return Amn;
}


//______________________________________________________________________________
double FourierBesselDecomp2D::ComputeImAmn(int m,int n,TH2D* h2d)
{
  // compute Im{Amn} coeficients
  
  double Amn = 0.;
  int nix = h2d->GetNbinsX();
  int niy = h2d->GetNbinsY();
  double xcm = 0.;//h2d->GetMean(1);
  double ycm = 0.;//h2d->GetMean(2);
  double dx = h2d->GetXaxis()->GetBinWidth(1);
  double dy = h2d->GetYaxis()->GetBinWidth(1);
  for(int ix=0;ix<nix;++ix) {
    for(int iy=0;iy<niy;++iy) {
      double x = h2d->GetXaxis()->GetBinCenter(ix+1) - xcm;
      double y = h2d->GetYaxis()->GetBinCenter(iy+1) - ycm;
      double fxy = h2d->GetBinContent(ix+1,iy+1);
      Amn += (1./(TMath::Pi()*r0*r0))*(fxy*Im_phi_mn(m,n,x,y)*dx*dy);
    }
  }
  
  return Amn;
}


//______________________________________________________________________________
double FourierBesselDecomp2D::func(double* x,double* Amn)
{
  // f(x,y) = Sum_mn{ A_mn*phi_mn(x,y) }
  
  double fxy = 0.;
  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  int maxdim = mdim*ndim;
  for(int m=m_min;m<=m_max;++m) {
    for(int n=n_min;n<=n_max;++n) {
      int indexRe = ndim*(m-m_min)+(n-n_min);
      int indexIm = ndim*(m-m_min)+(n-n_min) + maxdim;
      fxy += Amn[indexRe]*Re_phi_mn(m,n,x[0],x[1]) + Amn[indexIm]*Im_phi_mn(m,n,x[0],x[1]);
    }
  }

  return fxy;
}


//______________________________________________________________________________
TGraph2D* FourierBesselDecomp2D::ConvertToGraph2D(TH2D* h2d)
{
  // convert 2d histo to 2d graph
  
//  double xcm = h2d->GetMean(1);
//  double ycm = h2d->GetMean(2);
//  double deltax = h2d->GetXaxis()->GetBinWidth(1);
//  double deltay = h2d->GetYaxis()->GetBinWidth(1);

  TGraph2D* gr = new TGraph2D();
  int nx = h2d->GetNbinsX();
  int ny = h2d->GetNbinsY();
  for(int ix=0;ix<nx;++ix) {
    for(int iy=0;iy<ny;++iy) {
      int point = nx*ix + iy;
      double x = h2d->GetXaxis()->GetBinCenter(ix+1);
      double y = h2d->GetYaxis()->GetBinCenter(iy+1);
      double fxy = h2d->GetBinContent(ix+1,iy+1);
      cout << "point: " << point << "\tx: " << x << "\ty: " << y << "\tfxy: " << fxy << endl;
      gr->SetPoint(point,x,y,fxy);
    }
  }
  
  return gr;
}


//______________________________________________________________________________
void FourierBesselDecomp2D::DoDecomposition(TH2D* h2d)
{
  // calculate Fourier 2D decomposition

  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  int maxdim = mdim*ndim;

  if(hZeros) { delete hZeros; hZeros = 0; }
  hZeros = new TH2D("hZeros","#lambda_{m,n} vs m,n;m (angular);n (radial);#lambda_{m,n}",mdim,m_min-0.5,m_max+0.5,ndim,n_min-0.5,n_max+0.5);
  if(hReCoefs) { delete hReCoefs; hReCoefs = 0; }
  hReCoefs = new TH2D("hReCoefs","Re[A_{m,n}] vs m,n;m (angular);n (radial);Re[A_{m,n}]",mdim,m_min-0.5,m_max+0.5,ndim,n_min-0.5,n_max+0.5);
  if(hImCoefs) { delete hImCoefs; hImCoefs = 0; }
  hImCoefs = new TH2D("hImCoefs","Im[A_{m,n}] vs m,n;m (angular);n (radial);Im[A_{m,n}]",mdim,m_min-0.5,m_max+0.5,ndim,n_min-0.5,n_max+0.5);
  if(hCoefsNorm) { delete hCoefsNorm; hCoefsNorm = 0; }
  hCoefsNorm = new TH2D("hCoefsNorm","|A_{m,n}| vs m,n;m (angular);n (radial);|A_{m,n}|",mdim,m_min-0.5,m_max+0.5,ndim,n_min-0.5,n_max+0.5);

//  gr2d = (TGraph2D*)this->ConvertToGraph2D(h2d);
  
  int npar = 2.*maxdim;
  if(f2d) { delete f2d; f2d = 0; }
  f2d = new TF2("f2d",this,&FourierBesselDecomp2D::func,-8.,8.,-8.,8.,npar,"FourierBesselDecomp2D","func");
  //cout << endl;
  for(int m=m_min;m<=m_max;++m) {
    for(int n=n_min;n<=n_max;++n) {
      int indexRe = ndim*(m-m_min)+(n-n_min);
      int indexIm = ndim*(m-m_min)+(n-n_min) + maxdim;
      double ReAmn = ComputeReAmn(m,n,h2d);
      double ImAmn = ComputeImAmn(m,n,h2d);
//      cout << "im: " << im << "\tin: " << in << "\tindexRe: " << indexRe << "\tindexIm: " << indexIm << endl;
      f2d->SetParameter(indexRe,ReAmn); f2d->SetParName(indexRe,Form("Re[A_{%d,%d}]",m,n));
      f2d->SetParameter(indexIm,ImAmn); f2d->SetParName(indexIm,Form("Im[A_{%d,%d}]",m,n));
      //
      // fill histograms
      hReCoefs->Fill(m,n,ReAmn);
      hImCoefs->Fill(m,n,ImAmn);
      hCoefsNorm->Fill(m,n,TMath::Sqrt(ReAmn*ReAmn+ImAmn*ImAmn));
    }
  }
 
  // set npoints in TF2 to match the number of bins of h2d
  int nbx = h2d->GetNbinsX();
  int nby = h2d->GetNbinsY();
  f2d->SetNpx(nbx);
  f2d->SetNpy(nby);
 
}


//______________________________________________________________________________
double FourierBesselDecomp2D::GetL2()
{
  // return the L2 norm -- measure of the total mass of f (see arXiv:1204.5774v1)
  
  int npar = f2d->GetNpar();
  double L2 = 0.;
  for(int i=0;i<(npar/2);++i) {
    double ReAmn = f2d->GetParameter(i);
    double ImAmn = f2d->GetParameter(i+npar/2);
    L2 += pow(ReAmn,2.) + pow(ImAmn,2.);
  }
  
  return TMath::Sqrt(L2);
  
}

  
//______________________________________________________________________________
double FourierBesselDecomp2D::GetH1()
{
  // return the H1 norm -- measure of the wobblyness across the disk (see arXiv:1204.5774v1)
  
  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  int maxdim = mdim*ndim;
  double H1 = 0.;
  double l = 1.; // characteristic length scale -- see arXiv:1204.5774v1
  for(int m=m_min;m<=m_max;++m) {
    for(int n=n_min;n<=n_max;++n) {
      int indexRe = ndim*(m-m_min)+(n-n_min);
      int indexIm = ndim*(m-m_min)+(n-n_min) + maxdim;
      double ReAmn = f2d->GetParameter(indexRe);
      double ImAmn = f2d->GetParameter(indexIm);
      double lamb_mn = BesselZero(m,n);
      H1 += ((l*l)*(lamb_mn*lamb_mn) + 1)*(pow(ReAmn,2.) + pow(ImAmn,2.));
    }
  }
  
  return TMath::Sqrt(H1);
  
}

  
//______________________________________________________________________________
double FourierBesselDecomp2D::GetM1()
{
  // return the M1 norm -- measure of the angular gradient of f (see arXiv:1204.5774v1)
  
  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  int maxdim = mdim*ndim;
  double M1 = 0.;
  for(int m=m_min;m<=m_max;++m) {
    for(int n=n_min;n<=n_max;++n) {
      int indexRe = ndim*(m-m_min)+(n-n_min);
      int indexIm = ndim*(m-m_min)+(n-n_min) + maxdim;
      double ReAmn = f2d->GetParameter(indexRe);
      double ImAmn = f2d->GetParameter(indexIm);
//      double lamb_mn = BesselZero(m,n);
      M1 += (m*m)*(pow(ReAmn,2.) + pow(ImAmn,2.));
    }
  }
  
  return TMath::Sqrt(M1);
  
}

  
  //______________________________________________________________________________
void FourierBesselDecomp2D::Draw(TString opt)
{
  // plot function obtained

  TH2D* hframe = (TH2D*)f2d->GetHistogram();
  hframe->SetStats(false);
  hframe->SetTitle("2D Fourier-Bessel decomposition");
  hframe->GetXaxis()->SetTitle("x [fm]");
  hframe->GetXaxis()->SetTitleFont(132);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetXaxis()->SetTitleOffset(1.25);
  hframe->GetXaxis()->CenterTitle(true);
  hframe->GetXaxis()->SetLabelFont(132);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetXaxis()->SetLabelOffset(0.015);
  hframe->GetXaxis()->SetNdivisions(310);
  hframe->GetYaxis()->SetTitle("y [fm]");
  hframe->GetYaxis()->SetTitleFont(132);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleOffset(1.35);
  hframe->GetYaxis()->CenterTitle(true);
  hframe->GetYaxis()->SetLabelFont(132);
  hframe->GetYaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelOffset(0.015);
  hframe->GetYaxis()->SetNdivisions(310);
  hframe->GetZaxis()->SetTitleFont(132);
  hframe->GetZaxis()->SetTitleSize(0.06);
  hframe->GetZaxis()->SetTitleOffset(1.35);
  hframe->GetZaxis()->CenterTitle(true);
  hframe->GetZaxis()->SetLabelFont(132);
  hframe->GetZaxis()->SetLabelSize(0.05);
  hframe->GetZaxis()->SetLabelOffset(0.015);
  hframe->GetZaxis()->SetNdivisions(310);
  f2d->Draw(Form("%s",opt.Data()));
  
}


//______________________________________________________________________________
void FourierBesselDecomp2D::DrawReCoefsHisto(TString opt)
{
  // plot hReCoefs

  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  hReCoefs->SetStats(false);
  hReCoefs->GetXaxis()->SetTitleFont(132);
  hReCoefs->GetXaxis()->SetTitleSize(0.06);
  hReCoefs->GetXaxis()->SetTitleOffset(1.25);
  hReCoefs->GetXaxis()->CenterTitle(true);
  hReCoefs->GetXaxis()->SetLabelFont(132);
  hReCoefs->GetXaxis()->SetLabelSize(0.05);
  hReCoefs->GetXaxis()->SetLabelOffset(0.015);
  hReCoefs->GetXaxis()->SetNdivisions(100+mdim/2);
  hReCoefs->GetYaxis()->SetTitleFont(132);
  hReCoefs->GetYaxis()->SetTitleSize(0.06);
  hReCoefs->GetYaxis()->SetTitleOffset(1.35);
  hReCoefs->GetYaxis()->CenterTitle(true);
  hReCoefs->GetYaxis()->SetLabelFont(132);
  hReCoefs->GetYaxis()->SetLabelSize(0.05);
  hReCoefs->GetYaxis()->SetLabelOffset(0.015);
  hReCoefs->GetYaxis()->SetNdivisions(100+ndim);
  hReCoefs->GetZaxis()->SetTitleFont(132);
  hReCoefs->GetZaxis()->SetTitleSize(0.06);
  hReCoefs->GetZaxis()->SetTitleOffset(1.1);
  hReCoefs->GetZaxis()->CenterTitle(true);
  hReCoefs->GetZaxis()->SetLabelFont(132);
  hReCoefs->GetZaxis()->SetLabelSize(0.05);
  hReCoefs->GetZaxis()->SetLabelOffset(0.015);
  hReCoefs->GetZaxis()->SetNdivisions(310);
  hReCoefs->Draw(Form("%s",opt.Data()));
  
}


//______________________________________________________________________________
void FourierBesselDecomp2D::DrawImCoefsHisto(TString opt)
{
  // plot hImCoefs

  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  hImCoefs->SetStats(false);
  hImCoefs->GetXaxis()->SetTitleFont(132);
  hImCoefs->GetXaxis()->SetTitleSize(0.06);
  hImCoefs->GetXaxis()->SetTitleOffset(1.25);
  hImCoefs->GetXaxis()->CenterTitle(true);
  hImCoefs->GetXaxis()->SetLabelFont(132);
  hImCoefs->GetXaxis()->SetLabelSize(0.05);
  hImCoefs->GetXaxis()->SetLabelOffset(0.015);
  hImCoefs->GetXaxis()->SetNdivisions(100+mdim/2);
  hImCoefs->GetYaxis()->SetTitleFont(132);
  hImCoefs->GetYaxis()->SetTitleSize(0.06);
  hImCoefs->GetYaxis()->SetTitleOffset(1.35);
  hImCoefs->GetYaxis()->CenterTitle(true);
  hImCoefs->GetYaxis()->SetLabelFont(132);
  hImCoefs->GetYaxis()->SetLabelSize(0.05);
  hImCoefs->GetYaxis()->SetLabelOffset(0.015);
  hImCoefs->GetYaxis()->SetNdivisions(100+ndim);
  hImCoefs->GetZaxis()->SetTitleFont(132);
  hImCoefs->GetZaxis()->SetTitleSize(0.06);
  hImCoefs->GetZaxis()->SetTitleOffset(1.1);
  hImCoefs->GetZaxis()->CenterTitle(true);
  hImCoefs->GetZaxis()->SetLabelFont(132);
  hImCoefs->GetZaxis()->SetLabelSize(0.05);
  hImCoefs->GetZaxis()->SetLabelOffset(0.015);
  hImCoefs->GetZaxis()->SetNdivisions(310);
  hImCoefs->Draw(Form("%s",opt.Data()));
  
}


//______________________________________________________________________________
void FourierBesselDecomp2D::DrawCoefsNormHisto(TString opt)
{
  // plot hCoefsNorm

  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  hCoefsNorm->SetStats(false);
  hCoefsNorm->GetXaxis()->SetTitleFont(132);
  hCoefsNorm->GetXaxis()->SetTitleSize(0.06);
  hCoefsNorm->GetXaxis()->SetTitleOffset(1.25);
  hCoefsNorm->GetXaxis()->CenterTitle(true);
  hCoefsNorm->GetXaxis()->SetLabelFont(132);
  hCoefsNorm->GetXaxis()->SetLabelSize(0.05);
  hCoefsNorm->GetXaxis()->SetLabelOffset(0.015);
  hCoefsNorm->GetXaxis()->SetNdivisions(100+mdim/2);
  hCoefsNorm->GetYaxis()->SetTitleFont(132);
  hCoefsNorm->GetYaxis()->SetTitleSize(0.06);
  hCoefsNorm->GetYaxis()->SetTitleOffset(1.35);
  hCoefsNorm->GetYaxis()->CenterTitle(true);
  hCoefsNorm->GetYaxis()->SetLabelFont(132);
  hCoefsNorm->GetYaxis()->SetLabelSize(0.05);
  hCoefsNorm->GetYaxis()->SetLabelOffset(0.015);
  hCoefsNorm->GetYaxis()->SetNdivisions(100+ndim);
  hCoefsNorm->GetZaxis()->SetTitleFont(132);
  hCoefsNorm->GetZaxis()->SetTitleSize(0.06);
  hCoefsNorm->GetZaxis()->SetTitleOffset(1.1);
  hCoefsNorm->GetZaxis()->CenterTitle(true);
  hCoefsNorm->GetZaxis()->SetLabelFont(132);
  hCoefsNorm->GetZaxis()->SetLabelSize(0.05);
  hCoefsNorm->GetZaxis()->SetLabelOffset(0.015);
  hCoefsNorm->GetZaxis()->SetNdivisions(310);
  hCoefsNorm->Draw(Form("%s",opt.Data()));
  
}


//______________________________________________________________________________
void FourierBesselDecomp2D::DrawZerosHisto(TString opt)
{
  // plot hZeros

  int mdim = m_max-m_min+1;
  int ndim = n_max-n_min+1;
  hZeros->SetStats(false);
  hZeros->GetXaxis()->SetTitleFont(132);
  hZeros->GetXaxis()->SetTitleSize(0.06);
  hZeros->GetXaxis()->SetTitleOffset(1.25);
  hZeros->GetXaxis()->CenterTitle(true);
  hZeros->GetXaxis()->SetLabelFont(132);
  hZeros->GetXaxis()->SetLabelSize(0.05);
  hZeros->GetXaxis()->SetLabelOffset(0.015);
  hZeros->GetXaxis()->SetNdivisions(100+mdim);
  hZeros->GetYaxis()->SetTitleFont(132);
  hZeros->GetYaxis()->SetTitleSize(0.06);
  hZeros->GetYaxis()->SetTitleOffset(1.35);
  hZeros->GetYaxis()->CenterTitle(true);
  hZeros->GetYaxis()->SetLabelFont(132);
  hZeros->GetYaxis()->SetLabelSize(0.05);
  hZeros->GetYaxis()->SetLabelOffset(0.015);
  hZeros->GetYaxis()->SetNdivisions(100+ndim);
  hZeros->GetZaxis()->SetTitleFont(132);
  hZeros->GetZaxis()->SetTitleSize(0.06);
  hZeros->GetZaxis()->SetTitleOffset(1.1);
  hZeros->GetZaxis()->CenterTitle(true);
  hZeros->GetZaxis()->SetLabelFont(132);
  hZeros->GetZaxis()->SetLabelSize(0.05);
  hZeros->GetZaxis()->SetLabelOffset(0.015);
  hZeros->GetZaxis()->SetNdivisions(310);
  hZeros->Draw(Form("%s",opt.Data()));
  
}


//______________________________________________________________________________
void FourierBesselDecomp2D::Plot_phi_mn(int m,int n,double xmin,double xmax,double ymin,double ymax,TString comp)
{
  // plot phi_mn(x,y)

  TF2* fphi = 0;
  if(comp.Contains("Re")) fphi = new TF2("fphi","((ROOT::Math::cyl_bessel_j([0],TMath::Sqrt(x*x+y*y)*[1]))/(ROOT::Math::cyl_bessel_j(TMath::Abs([0])+1.,[2]*[1])))*(TMath::Cos([0]*TMath::ATan2(y,x)))",xmin,xmax,ymin,ymax);
  if(comp.Contains("Im")) fphi = new TF2("fphi","((ROOT::Math::cyl_bessel_j([0],TMath::Sqrt(x*x+y*y)*[1]))/(ROOT::Math::cyl_bessel_j(TMath::Abs([0])+1.,[2]*[1])))*(TMath::Sin([0]*TMath::ATan2(y,x)))",xmin,xmax,ymin,ymax);
  double lamb_mn = lambda[m-m_min][n-n_min]; // from the table BesselJmZeros.h
  fphi->SetParameter(0,double(m));
  fphi->SetParameter(1,lamb_mn);
  fphi->SetParameter(2,r0);
  TH2D* hframe = (TH2D*)fphi->GetHistogram();
  hframe->SetStats(false);
  hframe->SetTitle(Form("%s[#phi_{%d,%d}](x,y)",comp.Data(),m,n));
  hframe->GetXaxis()->SetTitle("x [fm]");
  hframe->GetXaxis()->SetTitleFont(132);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetXaxis()->SetTitleOffset(1.25);
  hframe->GetXaxis()->CenterTitle(true);
  hframe->GetXaxis()->SetLabelFont(132);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetXaxis()->SetLabelOffset(0.015);
  hframe->GetXaxis()->SetNdivisions(310);
  hframe->GetYaxis()->SetTitle("y [fm]");
  hframe->GetYaxis()->SetTitleFont(132);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleOffset(1.35);
  hframe->GetYaxis()->CenterTitle(true);
  hframe->GetYaxis()->SetLabelFont(132);
  hframe->GetYaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelOffset(0.015);
  hframe->GetYaxis()->SetNdivisions(310);
  hframe->GetZaxis()->SetTitleFont(132);
  hframe->GetZaxis()->SetTitleSize(0.06);
  hframe->GetZaxis()->SetTitleOffset(1.35);
  hframe->GetZaxis()->CenterTitle(true);
  hframe->GetZaxis()->SetLabelFont(132);
  hframe->GetZaxis()->SetLabelSize(0.05);
  hframe->GetZaxis()->SetLabelOffset(0.015);
  hframe->GetZaxis()->SetNdivisions(310);
  fphi->Draw("surf1");
   
}
