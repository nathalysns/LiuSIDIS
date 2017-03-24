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

double Esintheta_A(const double * x, const double * par);

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
  h0->GetXaxis()->SetTitle("sin#it{#theta_{#mu}}");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleFont(22);
  h0->GetXaxis()->SetTitleSize(0.055);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetTitle("sin#it{#theta_{A}}");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleFont(22);
  h0->GetYaxis()->SetTitleSize(0.055);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetRangeUser(-1.0, 1.0);

  TCanvas * c0 = new TCanvas("c0", "", 800, 800);
  c0->SetBottomMargin(0.15);
  c0->SetLeftMargin(0.15);

  if (opt == 1){//susy with d-EDM
    c0->Clear();
    double x[201];
    for (int i = 0; i < 201; i++){
      x[i] = -1.0 + 0.01 * i;
    }
    double md = 4.7e-3;//d quark mass in GeV
    double alphaS = 0.41 * M_PI;
    double Lambda = 1.0e3;//
    double tanbeta = 3.0;
    //par: dd, md, alpha_s, Lambda, tanbeta, 
    double par1[5] = {1.02e-24, md, alphaS, Lambda, tanbeta};//current gT + current EDMs
    double y1[201], e1[201];
    for (int i = 0; i < 201; i++){
      y1[i] = x[i] * par1[4];
      e1[i] = Esintheta_A(x, par1);
    }
    cout << e1[1] << endl;
    TGraphErrors * g1 = new TGraphErrors(201, x, y1, 0, e1);
    g1->SetLineColor(4);
    g1->SetFillColor(4);
    double par2[5] = {5.18e-25, md, alphaS, Lambda, tanbeta};//current gT + current EDMs
    double y2[201], e2[201];
    for (int i = 0; i < 201; i++){
      y2[i] = x[i] * par2[4];
      e2[i] = Esintheta_A(x, par2);
    }
    cout << e2[1] << endl;
    TGraphErrors * g2 = new TGraphErrors(201, x, y2, 0, e2);
    g2->SetLineColor(2);
    g2->SetFillColor(2);
    double par3[5] = {4.39e-25, md, alphaS, Lambda, tanbeta};//current gT + current EDMs
    double y3[201], e3[201];
    for (int i = 0; i < 201; i++){
      y3[i] = x[i] * par3[4];
      e3[i] = Esintheta_A(x, par3);
    }
    cout << e3[1] << endl;
    TGraphErrors * g3 = new TGraphErrors(201, x, y3, 0, e3);
    g3->SetLineColor(1);
    g3->SetFillColor(1);
    TLegend * l0 = new TLegend(0.15, 0.8, 0.55, 0.9);
    l0->AddEntry(g1, "#font[22]{current g_{T} + current d_{n} + current d_{p}}", "f");
    l0->AddEntry(g2, "#font[22]{future  g_{T} + current d_{n} + current d_{p}}", "f");
    l0->AddEntry(g3, "#font[22]{future  g_{T} + future d_{n} + current d_{p}}", "f");
    
    h0->DrawClone("axis");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g3->DrawClone("3same");
    h0->DrawClone("axissame");
    l0->DrawClone("same");
    c0->Print("susyd.pdf");
  }

  return 0;
}

//par: dd, md, alpha_s, Lambda, beta, 

double Esintheta_A(const double * x, const double * par){//
  double dd = par[0];//d quark edm / e, in unit cm
  dd = dd * 1.0e13 / 0.197327;//change to GeV^-1
  double md = par[1];
  double alpha_s = par[2];
  double Lambda = par[3];
  double EsinA = dd * (18.0 * M_PI * Lambda * Lambda / ((1.0/3.0) * alpha_s * md));
  return EsinA;
}
