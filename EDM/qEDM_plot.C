#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"

using namespace std;

double uedm_p(const double * x, const double * par);
double uedm_n(const double * x, const double * par);
double uedm_p_error(const double * x, const double * par);
double uedm_n_error(const double * x, const double * par);

double uedm_limit(const double * x, const double * par);
double dedm_limit(const double * x, const double * par);


int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./qEDM_plot <opt>" << endl;
    return -1;
  }

  int opt = atoi(argv[1]);

  TH1D * h0 = new TH1D("h0", "", 1, -1.0, 1.0);
  h0->SetStats(0);
  h0->SetTitle("");
  h0->GetXaxis()->SetNdivisions(5, 5, 0);
  h0->GetYaxis()->SetNdivisions(5, 5, 0);
  h0->GetXaxis()->SetTitle("");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleFont(22);
  h0->GetXaxis()->SetTitleSize(0.055);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetTitle("");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleFont(22);
  h0->GetYaxis()->SetTitleSize(0.055);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetLabelOffset(0.01);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetBottomMargin(0.15);
  c0->SetLeftMargin(0.15);

  if (opt == 1){//pEDM plot
    c0->Clear();
    double scale = 1.0e-24;
    h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-24} #it{e}#upointcm");
    h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-24} #it{e}#upointcm");
    h0->GetYaxis()->SetRangeUser(-1.2, 1.4);
    h0->GetYaxis()->SetNdivisions(6, 5, 0);
    double dd[400];
    for (int i = 0; i < 400; i++){
      dd[i] = -2.0 + 0.01 * i;
    }
    double par1[7] = {2.0e-25 / scale, 0.0, 0.413, 0.133, -0.229, 0.094, 0.002};//current pEDM + current gT
    double du1[400], eu1[400];
    for (int i = 0; i < 400; i++){
      du1[i] = uedm_p(&dd[i], par1);
      eu1[i] = uedm_p_error(&dd[i], par1);
    }
    TGraphErrors * g1 = new TGraphErrors(400, dd, du1, 0, eu1);
    g1->SetLineColor(4);
    g1->SetFillColor(4);
    double par2[7] = {2.0e-25 / scale, 0.0, 0.413, 0.018, -0.229, 0.008, -2.65e-5};//current pEDM + future gT
    double du2[400], eu2[400];
    for (int i = 0; i < 400; i++){
      du2[i] = uedm_p(&dd[i], par2);
      eu2[i] = uedm_p_error(&dd[i], par2);
    }
    TGraphErrors * g2 = new TGraphErrors(400, dd, du2, 0, eu2);
    g2->SetLineColor(2);
    g2->SetFillColor(2);
    double par3[7] = {2.0e-29 / scale, 0.0, 0.413, 0.018, -0.229, 0.008, -2.65e-5};//current pEDM + future gT
    double du3[400], eu3[400];
    for (int i = 0; i < 400; i++){
      du3[i] = uedm_p(&dd[i], par3);
      eu3[i] = uedm_p_error(&dd[i], par3);
    }
    TGraphErrors * g3 = new TGraphErrors(400, dd, du3, 0, eu3);
    g3->SetLineColor(1);
    g3->SetFillColor(1);
    TLegend * l0 = new TLegend(0.15, 0.78, 0.45, 0.9);
    l0->AddEntry(g1, "#font[22]{current g_{T} + current d_{p}}", "f");
    l0->AddEntry(g2, "#font[22]{future  g_{T} + current d_{p}}", "f");
    l0->AddEntry(g3, "#font[22]{future  g_{T} + future  d_{p}}", "f");
    
    h0->DrawClone("axis");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g3->DrawClone("3same");
    h0->DrawClone("axissame");
    l0->DrawClone("same");
    c0->Print("pedm.pdf");
  }

  if (opt == 2){//nEDM plot
    c0->Clear();
    double scale = 1.0e-25;
    h0->GetXaxis()->SetTitle("#it{d_{d}} / 10^{-25} #it{e}#upointcm");
    h0->GetYaxis()->SetTitle("#it{d_{u}} / 10^{-25} #it{e}#upointcm");
    h0->GetYaxis()->SetRangeUser(-3.0, 3.5);
    h0->GetYaxis()->SetNdivisions(7, 5, 0);
    double dd[400];
    for (int i = 0; i < 400; i++){
      dd[i] = -2.0 + 0.01 * i;
    }
    double par1[7] = {0.0, 2.1e-26 / scale, 0.413, 0.133, -0.229, 0.094, 0.002};//current nEDM + current gT
    double du1[400], eu1[400];
    for (int i = 0; i < 400; i++){
      du1[i] = uedm_n(&dd[i], par1);
      eu1[i] = uedm_n_error(&dd[i], par1);
    }
    TGraphErrors * g1 = new TGraphErrors(400, dd, du1, 0, eu1);
    g1->SetLineColor(4);
    g1->SetFillColor(4);
    double par2[7] = {0.0, 2.1e-26 / scale, 0.413, 0.018, -0.229, 0.008, -2.65e-5};//current nEDM + future gT
    double du2[400], eu2[400];
    for (int i = 0; i < 400; i++){
      du2[i] = uedm_n(&dd[i], par2);
      eu2[i] = uedm_n_error(&dd[i], par2);
    }
    TGraphErrors * g2 = new TGraphErrors(400, dd, du2, 0, eu2);
    g2->SetLineColor(2);
    g2->SetFillColor(2);
    double par3[7] = {0.0, 2.1e-28 / scale, 0.413, 0.018, -0.229, 0.008, -2.65e-5};//current nEDM + future gT
    double du3[400], eu3[400];
    for (int i = 0; i < 400; i++){
      du3[i] = uedm_n(&dd[i], par3);
      eu3[i] = uedm_n_error(&dd[i], par3);
    }
    TGraphErrors * g3 = new TGraphErrors(400, dd, du3, 0, eu3);
    g3->SetLineColor(1);
    g3->SetFillColor(1);

    TLegend * l0 = new TLegend(0.15, 0.78, 0.45, 0.9);
    l0->AddEntry(g1, "#font[22]{current g_{T} + current d_{n}}", "f");
    l0->AddEntry(g2, "#font[22]{future  g_{T} + current d_{n}}", "f");
    l0->AddEntry(g3, "#font[22]{future  g_{T} + future d_{n}}", "f");

    h0->DrawClone("axis");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g3->DrawClone("3same");
    h0->DrawClone("axissame");
    l0->DrawClone("same");
    c0->Print("nedm.pdf");
  }

  if (opt == 3){// qEDM limits
    double * xx;
    double par1[7] = {2.0e-25, 2.1e-26, 0.413, 0.133, -0.229, 0.094, 0.002};
    cout << "1: " << uedm_limit(xx, par1) << "  " << dedm_limit(xx, par1) << endl;
    double par2[7] = {2.0e-25, 2.1e-26, 0.413, 0.018, -0.229, 0.008, -2.65e-5};
    cout << "2: " << uedm_limit(xx, par2) << "  " << dedm_limit(xx, par2) << endl;
    double par3[7] = {2.0e-25, 2.1e-28, 0.413, 0.018, -0.229, 0.008, -2.65e-5};
    cout << "3: " << uedm_limit(xx, par3) << "  " << dedm_limit(xx, par3) << endl;
    double par4[7] = {2.0e-29, 2.1e-26, 0.413, 0.018, -0.229, 0.008, -2.65e-5};
    cout << "4: " << uedm_limit(xx, par4) << "  " << dedm_limit(xx, par4) << endl;
    double par5[7] = {2.0e-29, 2.1e-28, 0.413, 0.018, -0.229, 0.008, -2.65e-5};
    cout << "5: " << uedm_limit(xx, par5) << "  " << dedm_limit(xx, par5) << endl;
  }

  return 0;
}

//par: pedm, nedm, gT(u), eT(u), gT(d), eT(d), eT(ud)

double uedm_p(const double * x, const double * par){//
  double dd = x[0];//d quark edm
  double tu = par[2];
  double td = par[4];
  double du = (0.0 - td * dd) / tu;
  return du;
}

double uedm_n(const double * x, const double * par){//
  double dd = x[0];//d quark edm
  double tu = par[4];
  double td = par[2];
  double du = (0.0 - td * dd) / tu;
  return du;
}

double uedm_p_error(const double * x, const double * par){//
  double dd = x[0];//d quark edm
  double dp = par[0];//pedm limit
  double tu = par[2];
  double eu = par[3];
  double td = par[4];
  double ed = par[5];
  double eud = par[6];
  double edu2 =
    pow( (abs(dp) + abs(td * dd)) / tu, 2) * pow(eu / tu, 2)
    + pow( dd / tu, 2) * pow(ed, 2)
    + 2.0 * abs(dd * (dp + abs(td * dd)) / pow(tu, 3)) * eud;
  return sqrt(edu2) + abs(dp / tu);
}

double uedm_n_error(const double * x, const double * par){//
  double dd = x[0];//d quark edm
  double dn = par[1];//pedm limit
  double tu = par[4];
  double eu = par[5];
  double td = par[2];
  double ed = par[3];
  double eud = par[6];
  double edu2 =
    pow( (abs(dn) + abs(td * dd)) / tu, 2) * pow(eu / tu, 2)
    + pow( dd / tu, 2) * pow(ed, 2)
    + 2.0 * abs(dd * (dn + abs(td * dd)) / pow(tu, 3)) * eud;
  return sqrt(edu2) + abs(dn / tu);
}

double uedm_limit(const double * x, const double * par){//
  double dp = par[0];
  double dn = -par[1];
  double tu = par[2];
  double eu = par[3];
  double td = par[4];
  double ed = par[5];
  double eud = par[6];
  double Dt2 = tu * tu - td * td;
  double edu2 =
    pow( (-2.0 * tu * (tu * dp - td * dn)) / (Dt2 * Dt2) + dp / Dt2, 2) * pow(eu, 2)
    + pow( (2.0 * td * (tu * dp - td * dn)) / (Dt2 * Dt2) - dn / Dt2, 2) * pow(ed, 2)
    + 2.0 * ((-2.0 * tu * (tu * dp - td * dn)) / (Dt2 * Dt2) + dp / Dt2) * ((2.0 * td * (tu * dp - td * dn)) / (Dt2 * Dt2) - dn / Dt2) * eud;
  return sqrt(edu2) + abs( (tu * dp - td * dn) / Dt2);
}
 
double dedm_limit(const double * x, const double * par){//
  double dp = par[0];
  double dn = -par[1];
  double tu = par[2];
  double eu = par[3];
  double td = par[4];
  double ed = par[5];
  double eud = par[6];
  double Dt2 = tu * tu - td * td;
  double edu2 =
    pow( (-2.0 * tu * (tu * dn - td * dp)) / (Dt2 * Dt2) + dn / Dt2, 2) * pow(eu, 2)
    + pow( (2.0 * td * (tu * dn - td * dp)) / (Dt2 * Dt2) - dp / Dt2, 2) * pow(ed, 2)
    + 2.0 * ((-2.0 * tu * (tu * dn - td * dp)) / (Dt2 * Dt2) + dn / Dt2) * ((2.0 * td * (tu * dn - td * dp)) / (Dt2 * Dt2) - dp / Dt2) * eud;
  return sqrt(edu2) + abs( (tu * dn - td * dp) / Dt2);
}
