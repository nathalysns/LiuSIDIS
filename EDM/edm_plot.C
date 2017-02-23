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
double (*efx)(const double *, const double *);
double efx1(const double * x, const double * par);
double efx2(const double * x, const double * par);

int main(int argc, char * argv[]){

  if (argc < 2){
    printf("./edm_plot <opt>\n");
    return 1;
  }

  const double Xscale = 1e-25;
  const double Yscale = 1e-25;
  const double Xmin = -1.01;
  const double Xmax = 1.01;
  const double Ymin = -2.6;
  const double Ymax = 3.6;
  //base (blank)
  TH1D * h0 = new TH1D("h0", "", 1, Xmin, Xmax);
  h0->SetStats(0);
  h0->SetTitle("");
  h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-25} #it{e}#upointcm");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleFont(22);
  h0->GetXaxis()->SetTitleSize(0.055);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-25} #it{e}#upointcm");
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

  double x[100], y[100], err[100], err1[100], err2[100];
  double y1[100], y2[100];
  for (int i = 0; i < 100; i++)
    x[i] = (Xmin + (Xmax - Xmin) / 99.0 * i) * 1.0;
  
  const int opt = atoi(argv[1]);
  const double nedm68 = 1.82e-26 / Xscale;
  const double nedm90 = 2.99e-26 / Xscale;
  const double nedm95 = 3.57e-26 / Xscale;
  
  efx = efx1;
  if (argc > 2){
    efx = efx2;
  }
  
  if (opt == 1){//for SoLID
    c0->Clear();
    h0->DrawClone("AXIS");
    //par: ds, dn, tu, td, ts, edn, etu, etd, ets
    //KPSY15
    double par1[9] = {0.0, 0.0, 0.413, -0.229, 0.0, nedm90, 0.133, 0.094, 0.0};
    for (int i = 0; i < 100; i++){
      y[i] = fx(&x[i], par1);
      err[i] = efx(&x[i], par1);
    }
    TGraphErrors * g1 = new TGraphErrors(100, x, y, 0, err);
    g1->SetLineColorAlpha(4, 1.0);
    g1->SetFillColorAlpha(4, 1.0);
    //SoLID
    double par2[9] = {0.0, 0.0, 0.413, -0.229, 0.0, nedm90, 0.018, 0.008, 0.0};//SoLID
    for (int i = 0; i < 100; i++){
      y[i] = fx(&x[i], par2);
      err[i] = efx(&x[i], par2);
    }
    TGraphErrors * g2 = new TGraphErrors(100, x, y, 0, err);
    g2->SetLineColorAlpha(2, 1.0);
    g2->SetFillColorAlpha(2, 1.0);
    //SoLID + nEDM
    double par3[9] = {0.0, 0.0, 0.413, -0.229, 0.0, 0.01 * nedm90, 0.018, 0.008, 0.0};//SoLID
    for (int i = 0; i < 100; i++){
      y[i] = fx(&x[i], par3);
      err[i] = efx(&x[i], par3);
    }
    TGraphErrors * g3 = new TGraphErrors(100, x, y, 0, err);
    g3->SetLineColorAlpha(1, 1.0);
    g3->SetFillColorAlpha(1, 1.0);
    TLegend * l0 = new TLegend(0.15, 0.72, 0.45, 0.9);
    l0->AddEntry(g1, "#font[22]{Kang et al. (2016)}", "f");
    l0->AddEntry(g2, "#font[22]{SoLID projection}", "f");
    l0->AddEntry(g3, "#font[22]{SoLID + future nEDM}", "f");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g3->DrawClone("3same");
    h0->DrawClone("axissame");
    l0->Draw("same");
    c0->Print("solid_edm.pdf");
    
    for (int i = 0; i < 100; i++){
      y[i] = -x[i];
    }
    TGraph * gl = new TGraph(100, x, y);
    gl->SetLineColorAlpha(1, 0.8);
    gl->SetLineWidth(2);
    gl->SetLineStyle(2);
    gl->DrawClone("Lsame");
    //c0->Print("solid_edm_1.pdf");
    double par4[9] = {0.0, 0.0, 0.413, -0.229, 0.0, nedm90, 0.0, 0.0, 0.0};//SoLID
    for (int i = 0; i < 100; i++){
      y[i] = fx(&x[i], par4);
      err[i] = efx(&x[i], par4);
      err1[i] = y[i] + err[i];
      err2[i] = y[i] - err[i];
    }
    TGraph * limit1 = new TGraph(100, x, err1);
    limit1->SetLineColor(5);
    limit1->SetLineWidth(2);
    limit1->SetLineStyle(5);
    TGraph * limit2 = new TGraph(100, x, err2);
    limit2->SetLineColor(5);
    limit2->SetLineWidth(2);
    limit2->SetLineStyle(5);
    limit1->DrawClone("Lsame");
    limit2->DrawClone("Lsame");
    c0->Print("solid_edm_1.pdf");
  }

  if (opt == 2){//lattice
    c0->Clear();
    h0->DrawClone("AXIS");
    //par: ds, dn, tu, td, ts, edn, etu, etd, ets
    //PNDME16 PRD94(2016)054508
    double par1[9] = {0.0, 0.0, 0.792, -0.194, 0.0, nedm68, 0.042, 0.014, 0.0};
    for (int i = 0; i < 100; i++){
      y[i] = fx(&x[i], par1);
      err[i] = efx(&x[i], par1);
      y1[i] = y[i] + err[i];
      y2[i] = y[i] - err[i];
    }
    TGraphErrors * g1 = new TGraphErrors(100, x, y, 0, err);
    g1->SetLineColorAlpha(1, 1.0);
    g1->SetFillColorAlpha(1, 1.0);
    TGraph * g1a = new TGraph(100, x, y1);
    g1a->SetLineColorAlpha(1, 1.0);
    g1a->SetLineWidth(2);
    g1a->SetLineStyle(2);
    TGraph * g1b = new TGraph(100, x, y2);
    g1b->SetLineColorAlpha(1, 1.0);
    g1b->SetLineWidth(2);
    g1b->SetLineStyle(2);
    //ETMC15 PRD92(2015)114513
    double par2[9] = {0.0, 0.0, 0.791, -0.236, 0.0, nedm68, 0.053, 0.033, 0.0};//SoLID
    for (int i = 0; i < 100; i++){
      y[i] = fx(&x[i], par2);
      err[i] = efx(&x[i], par2);
    }
    TGraphErrors * g2 = new TGraphErrors(100, x, y, 0, err);
    g2->SetLineColorAlpha(2, 1.0);
    g2->SetFillColorAlpha(2, 1.0);
    //DSE PRD91(2015)074004
    double par3[9] = {0.0, 0.0, 0.55, -0.11, 0.0, nedm68, 0.08, 0.02, 0.0};//SoLID
    for (int i = 0; i < 100; i++){
      y[i] = fx(&x[i], par3);
      err[i] = efx(&x[i], par3);
    }
    TGraphErrors * g3 = new TGraphErrors(100, x, y, 0, err);
    g3->SetLineColorAlpha(4, 1.0);
    g3->SetFillColorAlpha(4, 1.0);
    //g3->DrawClone("3same");
    TLegend * l0 = new TLegend(0.15, 0.72, 0.45, 0.9);
    l0->AddEntry(g1, "#font[22]{LQCD Bhattacharya et al. (2016)}", "f");
    l0->AddEntry(g2, "#font[22]{LQCD Abdel-Rehim et al. (2015)}", "f");
    l0->AddEntry(g3, "#font[22]{DSE Pitschmann et al. (2015)}", "f");
    g3->DrawClone("3same");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g1a->DrawClone("Lsame");
    g1b->DrawClone("Lsame");
    h0->DrawClone("axissame");
    l0->DrawClone("same");
    c0->Print("lattice_dse_edm.pdf");
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
  double du = fx(x, par);
  double edu2 = 
    max(pow((edn - tu * dd - ts * ds) / td, 2), pow((-edn - tu * dd - ts * ds) / td, 2)) * pow(etd / td, 2)
    + pow( etu * dd / td, 2)
    + pow( ets * ds / td, 2);
  return sqrt(edu2) + sqrt(pow(edn / td, 2));
}

double efx2(const double * x, const double * par){
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
  double du = fx(x, par);
  double edu2 = 
    max(pow((dn - tu * dd - ts * ds) / td, 2), pow((dn - tu * dd - ts * ds) / td, 2)) * pow(etd / td, 2)
    + pow( etu * dd / td, 2)
    + pow( ets * ds / td, 2)
    + pow( edn / td, 2);
  return sqrt(edu2);
}
