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
  h0->GetXaxis()->SetTitle("sin#it{#theta_{A}}");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleFont(22);
  h0->GetXaxis()->SetTitleSize(0.055);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetTitle("sin#it{#theta_{#mu}}");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleFont(22);
  h0->GetYaxis()->SetTitleSize(0.055);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetLabelOffset(0.01);
  h0->GetYaxis()->SetRangeUser(-1.0, 1.0);

  TCanvas * c0 = new TCanvas("c0", "", 650, 750);
  c0->SetTopMargin(0.20);
  c0->SetBottomMargin(0.15);
  c0->SetLeftMargin(0.15);

  if (opt == 1){//susy with d-EDM
    c0->Clear();
    double x[201];
    for (int i = 0; i < 201; i++){
      x[i] = -1.0 + 0.01 * i;
    }
    double mu = 2.2e-3;//u quark mass in GeV
    double md = 4.7e-3;//d quark mass in GeV
    double alphaS = 0.41 * M_PI;
    double Lambda = 1.0e3;//
    const double unit = 1.0e13 / 0.197327;
    double du, dd;
    //par: dd, md, alpha_s, Lambda, tanbeta,
    double factor1 = 18.0 * M_PI * Lambda * Lambda / alphaS;
    du = 1.3e-24 * unit;
    dd = 1.0e-24 * unit;
    double y1[201], e1[201];
    for (int i = 0; i < 201; i++){
      y1[i] = x[i];
      e1[i] = sqrt((abs(x[i]) + abs(2.0 * dd / md) * factor1) * (abs(x[i]) + abs(3.0 * du / (2.0 * mu)) * factor1)) - abs(x[i]);
    }
    cout << e1[100] << endl;
    TGraphErrors * g1 = new TGraphErrors(201, x, y1, 0, e1);
    g1->SetLineColor(4);
    g1->SetFillColor(4);
    g1->SetMarkerStyle(8);
    g1->SetMarkerSize(1.1);
    g1->SetMarkerColor(4);
    //
    du = 8.0e-25 * unit;
    dd = 5.2e-25 * unit;
    double y2[201], e2[201];
    for (int i = 0; i < 201; i++){
      y2[i] = x[i];
      e2[i] = sqrt((abs(x[i]) + abs(2.0 * dd / md) * factor1) * (abs(x[i]) + abs(3.0 * du / (2.0 * mu)) * factor1)) - abs(x[i]);
    }
    cout << e2[100] << endl;
    TGraphErrors * g2 = new TGraphErrors(201, x, y2, 0, e2);
    g2->SetLineColor(2);
    g2->SetFillColor(2);
    g2->SetMarkerStyle(8);
    g2->SetMarkerSize(1.1);
    g2->SetMarkerColor(2);
    //
    du = 5.4e-28 * unit;
    dd = 8.4e-28 * unit;
    double y3[201], e3[201];
    for (int i = 0; i < 201; i++){
      y3[i] = x[i];
      e3[i] = sqrt((abs(x[i]) + abs(2.0 * dd / md) * factor1) * (abs(x[i]) + abs(3.0 * du / (2.0 * mu)) * factor1)) - abs(x[i]);
    }
    cout << e3[100] << endl;
    TGraphErrors * g3 = new TGraphErrors(201, x, y3, 0, e3);
    g3->SetLineColor(1);
    g3->SetFillColor(1);
    g3->SetMarkerStyle(8);
    g3->SetMarkerSize(1.1);
    g3->SetMarkerColor(1);
    factor1 = factor1 * 900.0;
    double y4[201], e4[201];
    double y4a[201], y4b[201];
    for (int i = 0; i < 201; i++){
      y4[i] = x[i];
      e4[i] = sqrt((abs(x[i]) + abs(2.0 * dd / md) * factor1) * (abs(x[i]) + abs(3.0 * du / (2.0 * mu)) * factor1)) - abs(x[i]);
      y4a[i] = y4[i] + e4[i];
      y4b[i] = y4[i] - e4[i];
    }
    TGraph * g4a = new TGraph(201, x, y4a);
    TGraph * g4b = new TGraph(201, x, y4b);
    g4a->SetLineColor(5);
    g4a->SetLineStyle(2);
    g4a->SetLineWidth(2);
    g4b->SetLineColor(5);
    g4b->SetLineStyle(2);
    g4b->SetLineWidth(2);

    TGraphErrors * g4 = new TGraphErrors(1);
    g4->SetLineColor(kYellow+1);
    g4->SetLineStyle(2);
    g4->SetLineWidth(2);
    g4->SetFillColor(0);
    g4->SetMarkerStyle(8);
    g4->SetMarkerSize(1.1);
    g4->SetMarkerColor(5);

    TLegend * l0 = new TLegend(0.15, 0.8, 0.9, 0.88);
    l0->SetNColumns(2);
    l0->AddEntry(g1, "#font[22]{current g_{T} + current d_{N} #Lambda=1TeV   }", "f");
    l0->AddEntry(g2, "#font[22]{future   g_{T} + current d_{N} #Lambda=1TeV   }", "f");
    l0->AddEntry(g3, "#font[22]{future   g_{T} + future   d_{N} #Lambda=1TeV   }", "f");
    l0->AddEntry(g4, "#font[22]{future   g_{T} + future   d_{N} #Lambda=30TeV   }", "f");

    h0->DrawClone("axis");
    g1->DrawClone("3same");
    g2->DrawClone("3same");
    g3->DrawClone("3lsame");
    g4a->DrawClone("lsame");
    g4b->DrawClone("lsame");
    h0->DrawClone("axissame");
    l0->DrawClone("same");
    c0->Print("susy.pdf");
  }

  return 0;
}
