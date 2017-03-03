#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

using namespace std;

double fx(const double * x, const double * par);
double fxp(const double * x, const double * par);
double (*efx)(const double *, const double *);
double efx1(const double * x, const double * par);
double efx2(const double * x, const double * par);

int main(int argc, char * argv[]){

  if (argc < 2){
    printf("./edm_plot <opt>\n");
    return 1;
  }

  const double Xmin = -1.01;
  const double Xmax = 1.01;
  const double Ymin = -2.6;
  const double Ymax = 3.6;
  //base (blank)
  TH1D * h0 = new TH1D("h0", "", 1, Xmin, Xmax);
  h0->SetStats(0);
  h0->SetTitle("");
  h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-23} #it{e}#upointcm");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleFont(22);
  h0->GetXaxis()->SetTitleSize(0.055);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-23} #it{e}#upointcm");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleFont(22);
  h0->GetYaxis()->SetTitleSize(0.055);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetRangeUser(Ymin, Ymax);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetBottomMargin(0.15);
  c0->SetLeftMargin(0.15);

  double x[1000], y[1000], err[1000], err1[1000], err2[1000];
  double y1[1000], y2[1000];
  for (int i = 0; i < 1000; i++)
    x[i] = (Xmin + (Xmax - Xmin) / 99.0 * i) * 1.0;
  
  const int opt = atoi(argv[1]);
  double nedm95 = 1.6e-26;
  double pedm95 = 2.0e-25;

  
  if (opt == 1){//for SoLID
    c0->Clear();
    nedm95 = 1.6e-26;
    double Xscale = 1.e-25;
    h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-25} #it{e}#upointcm");
    h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-25} #it{e}#upointcm");
    h0->DrawClone("AXIS");
    //par: ds, dn, tu, td, ts, edn, etu, etd, ets
    efx = efx1;//neutron
    //KPSY15
    double par1[10] = {0.0, 0.0, 0.413, -0.229, 0.0, nedm95/Xscale, 0.133, 0.094, 0.0, 0.002};
    for (int i = 0; i < 1000; i++){
      y[i] = fx(&x[i], par1);
      err[i] = efx(&x[i], par1);
    }
    TGraphErrors * g1 = new TGraphErrors(1000, x, y, 0, err);
    g1->SetLineColorAlpha(4, 1.0);
    g1->SetFillColorAlpha(4, 1.0);
    //SoLID
    double par2[10] = {0.0, 0.0, 0.413, -0.229, 0.0, nedm95/Xscale, 0.018, 0.008, 0.0, -2.65e-5};//SoLID
    for (int i = 0; i < 1000; i++){
      y[i] = fx(&x[i], par2);
      err[i] = efx(&x[i], par2);
    }
    TGraphErrors * g2 = new TGraphErrors(1000, x, y, 0, err);
    g2->SetLineColorAlpha(2, 1.0);
    g2->SetFillColorAlpha(2, 1.0);
    //SoLID + nEDM
    double par3[10] = {0.0, 0.0, 0.413, -0.229, 0.0, 0.01 * nedm95/Xscale, 0.018, 0.008, 0.0, -2.65e-5};//SoLID
    for (int i = 0; i < 1000; i++){
      y[i] = fx(&x[i], par3);
      err[i] = efx(&x[i], par3);
    }
    TGraphErrors * g3 = new TGraphErrors(1000, x, y, 0, err);
    g3->SetLineColorAlpha(1, 1.0);
    g3->SetFillColorAlpha(1, 1.0);
    TLegend * l0 = new TLegend(0.15, 0.78, 0.45, 0.9);
    l0->AddEntry(g1, "#font[22]{KPSY2016 + current nEDM}", "f");
    l0->AddEntry(g2, "#font[22]{SoLID + current nEDM}", "f");
    l0->AddEntry(g3, "#font[22]{SoLID + future nEDM}", "f");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g3->DrawClone("3same");
    h0->DrawClone("axissame");
    l0->DrawClone("same");
    c0->Print("nEDM3.pdf");
  }

  if (opt == 2){//for SoLID
    c0->Clear();
    nedm95 = 3.0e-26;
    double Xscale = 1.e-25;
    h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-25} #it{e}#upointcm");
    h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-25} #it{e}#upointcm");
    h0->DrawClone("AXIS");
    //par: ds, dn, tu, td, ts, edn, etu, etd, ets
    efx = efx1;//neutron
    //KPSY15
    double par1[10] = {0.0, 0.0, 0.413, -0.229, 0.0, nedm95/Xscale, 0.133, 0.094, 0.0, 0.002};
    for (int i = 0; i < 1000; i++){
      y[i] = fx(&x[i], par1);
      err[i] = efx(&x[i], par1);
    }
    TGraphErrors * g1 = new TGraphErrors(1000, x, y, 0, err);
    g1->SetLineColorAlpha(4, 1.0);
    g1->SetFillColorAlpha(4, 1.0);
    //SoLID
    double par2[10] = {0.0, 0.0, 0.413, -0.229, 0.0, nedm95/Xscale, 0.018, 0.008, 0.0, -2.65e-5};//SoLID
    for (int i = 0; i < 1000; i++){
      y[i] = fx(&x[i], par2);
      err[i] = efx(&x[i], par2);
    }
    TGraphErrors * g2 = new TGraphErrors(1000, x, y, 0, err);
    g2->SetLineColorAlpha(2, 1.0);
    g2->SetFillColorAlpha(2, 1.0);
    //SoLID + nEDM
    double par3[10] = {0.0, 0.0, 0.413, -0.229, 0.0, 0.01 * nedm95/Xscale, 0.018, 0.008, 0.0, -2.65e-5};//SoLID
    for (int i = 0; i < 1000; i++){
      y[i] = fx(&x[i], par3);
      err[i] = efx(&x[i], par3);
    }
    TGraphErrors * g3 = new TGraphErrors(1000, x, y, 0, err);
    g3->SetLineColorAlpha(1, 1.0);
    g3->SetFillColorAlpha(1, 1.0);
    TLegend * l0 = new TLegend(0.15, 0.78, 0.45, 0.9);
    l0->AddEntry(g1, "#font[22]{KPSY2016 + current nEDM}", "f");
    l0->AddEntry(g2, "#font[22]{SoLID + current nEDM}", "f");
    l0->AddEntry(g3, "#font[22]{SoLID + future nEDM}", "f");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g3->DrawClone("3same");
    h0->DrawClone("axissame");
    l0->DrawClone("same");
    c0->Print("nEDM4.pdf");
  }

  if (opt == 3){
    c0->Clear();
    pedm95 = 2.0e-25;
    double Xscale = 1.e-24;
    h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-24} #it{e}#upointcm");
    h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-24} #it{e}#upointcm");
    h0->DrawClone("AXIS");
    //par: ds, dn, tu, td, ts, edn, etu, etd, ets
    efx = efx2;
    //KPSY15 + pEDM
    double par4[10] = {0.0, 0.0, 0.413, -0.229, 0.0, pedm95/Xscale, 0.133, 0.094, 0.0, 0.002};
    for (int i = 0; i < 1000; i++){
      y[i] = fxp(&x[i], par4);
      err[i] = efx(&x[i], par4);
    }
    TGraphErrors * g4 = new TGraphErrors(1000, x, y, 0, err);
    g4->SetLineColorAlpha(4, 1.0);
    g4->SetFillColorAlpha(4, 1.0);
    //SoLID + pEDM
    double par5[10] = {0.0, 0.0, 0.413, -0.229, 0.0, pedm95/Xscale, 0.018, 0.008, 0.0, -2.65e-5};//SoLID
    for (int i = 0; i < 1000; i++){
      y[i] = fxp(&x[i], par5);
      err[i] = efx(&x[i], par5);
    }
    TGraphErrors * g5 = new TGraphErrors(1000, x, y, 0, err);
    g5->SetLineColorAlpha(2, 1.0);
    g5->SetFillColorAlpha(2, 1.0);
  
    TLegend * l1 = new TLegend(0.15, 0.78, 0.45, 0.9);
    l1->AddEntry(g4, "#font[22]{KPSY2016 + current pEDM}", "f");
    l1->AddEntry(g5, "#font[22]{SoLID + current pEDM}", "f");
    g4->DrawClone("3same");
    g5->DrawClone("3same");
    h0->DrawClone("axissame");
    l1->DrawClone("same");
    c0->Print("pEDM3.pdf");
  }

  if (opt == 4){
    c0->Clear();
    pedm95 = 5.4e-24;
    double Xscale = 1.e-23;
    h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-23} #it{e}#upointcm");
    h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-23} #it{e}#upointcm");
    h0->DrawClone("AXIS");
    //par: ds, dn, tu, td, ts, edn, etu, etd, ets
    efx = efx2;
    //KPSY15 + pEDM
    double par4[10] = {0.0, 0.0, 0.413, -0.229, 0.0, pedm95/Xscale, 0.133, 0.094, 0.0, 0.002};
    for (int i = 0; i < 1000; i++){
      y[i] = fxp(&x[i], par4);
      err[i] = efx(&x[i], par4);
    }
    TGraphErrors * g4 = new TGraphErrors(1000, x, y, 0, err);
    g4->SetLineColorAlpha(4, 1.0);
    g4->SetFillColorAlpha(4, 1.0);
    //SoLID + pEDM
    double par5[10] = {0.0, 0.0, 0.413, -0.229, 0.0, pedm95/Xscale, 0.018, 0.008, 0.0, -2.65e-5};//SoLID
    for (int i = 0; i < 1000; i++){
      y[i] = fxp(&x[i], par5);
      err[i] = efx(&x[i], par5);
    }
    TGraphErrors * g5 = new TGraphErrors(1000, x, y, 0, err);
    g5->SetLineColorAlpha(2, 1.0);
    g5->SetFillColorAlpha(2, 1.0);
  
    TLegend * l1 = new TLegend(0.15, 0.78, 0.45, 0.9);
    l1->AddEntry(g4, "#font[22]{KPSY2016 + current pEDM}", "f");
    l1->AddEntry(g5, "#font[22]{SoLID + current pEDM}", "f");
    g4->DrawClone("3same");
    g5->DrawClone("3same");
    h0->DrawClone("axissame");
    l1->DrawClone("same");
    c0->Print("pEDM4.pdf");
  }


  return 0;
}

double fx(const double * x, const double * par){
  double dd = x[0];//d quark edm
  double ds = par[0];//s quark edm
  double dn = par[1];//neutron edm
  double tu = par[2];//u quark tensor charge
  double td = par[3];//d quark tensor charge
  double ts = par[4];//s quark tensor charge
  double du = (dn - tu * dd - ts * ds) / td;
  return du;
}

double fxp(const double * x, const double * par){
  double dd = x[0];//d quark edm
  double ds = par[0];//s quark edm
  double dp = par[1];//proton edm
  double tu = par[2];//u quark tensor charge
  double td = par[3];//d quark tensor charge
  double ts = par[4];//s quark tensor charge
  double du = (dp - td * dd - ts * ds) / tu;
  return du;
}
  
double efx1(const double * x, const double * par){
  double dd = x[0];//d quark edm
  double ds = par[0];//s quark edm
  double dn = par[1];//neutron edm
  double tu = par[2];//u quark tensor charge
  double td = par[3];//d quark tensor charge
  double ts = par[4];//s quark tensor charge
  double edn = par[5];//error of dn
  double etu = par[6];//error of tu
  double etd = par[7];//error of td
  double ets = par[8];//error of ts
  double etud = par[9];//error correlation term tutd
  double du = fx(x, par);
  double edu2 = 
    max(pow((edn - tu * dd - ts * ds) / td, 2), pow((-edn - tu * dd - ts * ds) / td, 2)) * pow(etd / td, 2)
    + pow( etu * dd / td, 2)
    + pow( ets * ds / td, 2)
    + 2.0 * max((edn - tu * dd - ts * ds) * dd / pow(td, 3) * etud,
		(-edn - tu * dd - ts * ds) * dd / pow(td, 3) * etud);
  return sqrt(edu2) + sqrt(pow(edn / td, 2));
}

double efx2(const double * x, const double * par){
  double dd = x[0];//d quark edm
  double ds = par[0];//s quark edm
  double dp = par[1];//proton edm
  double tu = par[2];//u quark tensor charge
  double td = par[3];//d quark tensor charge
  double ts = par[4];//s quark tensor charge
  double edp = par[5];//error of dp
  double etu = par[6];//error of tu
  double etd = par[7];//error of td
  double ets = par[8];//error of ts
  double etud = par[9];//error correlation tutd
  double du = fxp(x, par);
  double edu2 = 
    max(pow((edp - td * dd - ts * ds) / tu, 2), pow((-edp - td * dd - ts * ds) / tu, 2)) * pow(etu / tu, 2)
    + pow( etd * dd / tu, 2)
    + pow( ets * ds / tu, 2)
    + 2.0 * max((edp - td * dd - ts * ds) * dd / pow(tu, 3) * etud,
		(-edp - td * dd - ts * ds) * dd / pow(tu, 3) * etud);
  return sqrt(edu2) + sqrt(pow(edp / tu, 2));
}
