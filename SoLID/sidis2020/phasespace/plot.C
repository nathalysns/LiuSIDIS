#include <iostream>
#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

int SetStyle(TH2D * h){
  h->SetStats(0);
  h->GetXaxis()->CenterTitle(true);
  h->GetXaxis()->SetTitle("x");
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetLabelSize(0.055);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetRangeUser(1e-3, 1.0);
  h->GetYaxis()->CenterTitle(true);
  h->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetLabelSize(0.055);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetRangeUser(0.5, 1000);
  return 0;
}
  
int main(const int argc, const char * argv[]){

  gStyle->SetOptStat(0);

  TFile * fs = new TFile("qt-3he.root", "r");
  TH2D * h0 = (TH2D *) fs->Get("11GeV_qt<Q");
  TH2D * h1 = (TH2D *) fs->Get("11GeV_qt<0.33Q");
  TH2D * h2 = (TH2D *) fs->Get("11GeV_qt<0.2Q");
  TH2D * g0 = (TH2D *) fs->Get("8.8GeV_qt<Q");
  TH2D * g1 = (TH2D *) fs->Get("8.8GeV_qt<0.33Q");
  TH2D * g2 = (TH2D *) fs->Get("8.8GeV_qt<0.2Q");
  TFile * feic = new TFile("qt-eic.root", "r");
  TH2D * e0 = (TH2D *) feic->Get("qt<Q");
  TH2D * e1 = (TH2D *) feic->Get("qt<0.33Q");
  TH2D * e2 = (TH2D *) feic->Get("qt<0.2Q");
  

  double m = h0->GetMaximum();
  h1->SetMaximum(m);
  h2->SetMaximum(m);

  SetStyle(h0);
  SetStyle(h1);
  SetStyle(h2);
  SetStyle(g0);
  SetStyle(g1);
  SetStyle(g2);
  SetStyle(e0);
  SetStyle(e1);
  SetStyle(e2);

  e0->SetMarkerColor(4);
  e1->SetMarkerColor(4);
  e2->SetMarkerColor(4);
  h0->SetMarkerColor(2);
  h1->SetMarkerColor(2);
  h2->SetMarkerColor(2);
  g0->SetMarkerColor(2);
  g1->SetMarkerColor(2);
  g2->SetMarkerColor(2);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetLeftMargin(0.15);
  c0->SetBottomMargin(0.15);
  c0->SetLogx();
  c0->SetLogy();
  c0->SetLogz();

  e0->DrawClone("axis");
  h0->DrawClone("same");
  g0->DrawClone("same");
  e0->DrawClone("same");

  c0->Print("qt1.pdf");

  e1->DrawClone("axis");
  h1->DrawClone("same");
  g1->DrawClone("same");
  e1->DrawClone("same");

  c0->Print("qt2.pdf");

  e2->DrawClone("axis");
  h2->DrawClone("same");
  g2->DrawClone("same");
  e2->DrawClone("same");

  c0->Print("qt3.pdf");


  return 0;
}
  
  
