void drawIniConProfile(TString funcName="gaus2d",double b = 2.0,double r=7.5) {
 
  gROOT->ProcessLine(Form(".L %s.cpp",funcName.Data()));
  gROOT->ProcessLine(".L loadStyle.C");
  loadStyle();

  TF2* fIniCon = new TF2(funcName.Data(),myinicon,-10.,10.,-10.,10.,5);
  fIniCon->SetTitle("");
  double par[] = {1.,0.,0.,(r-b/2.)/3.,r/3.};
  int npar = sizeof(par)/sizeof(double);
  for(int i=0;i<npar;++i) fIniCon->SetParameter(i,par[i]);

  TCanvas* cc = new TCanvas("cc","Initial Energy Density Transverse Profile",10,10,1200,600);
  cc->Divide(2,1,0.0001,0.0001);
  cc->cd(2);
  cc->cd(2)->SetGridx();
  cc->cd(2)->SetGridy();
  fIniCon->GetXaxis()->SetTitle("x [fm]");
  fIniCon->GetXaxis()->SetTitleOffset(2.0);
  fIniCon->GetXaxis()->CenterTitle(true);
  fIniCon->GetXaxis()->SetNdivisions(305);
  fIniCon->GetYaxis()->SetTitle("y [fm]");
  fIniCon->GetYaxis()->SetTitleOffset(2.0);
  fIniCon->GetYaxis()->CenterTitle(true);
  fIniCon->GetYaxis()->SetNdivisions(305);
  fIniCon->GetZaxis()->SetNdivisions(310);
  fIniCon->Draw("contz");

  cc->cd(1);
  fIniCon->Draw("surf1fb");

}

inline double myinicon(double *x,double *par) { return inicon(x,2,par); }
