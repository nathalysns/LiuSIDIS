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
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetLabelSize(0.055);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetRangeUser(0., 1.);
  h->GetYaxis()->CenterTitle(true);
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetLabelSize(0.055);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetRangeUser(0.0, 2.0);
  return 0;
}
  
int main(const int argc, const char * argv[]){

  gStyle->SetOptStat(0);

  //TFile * fs = new TFile("zpt-3he.root", "r");
  TFile * fs = new TFile("zpt-3he-all.root", "r");
  TH2D * h0 = (TH2D *) fs->Get("11GeV_All");
  TH2D * h1 = (TH2D *) fs->Get("11GeV_R<0.4");
  TH2D * h2 = (TH2D *) fs->Get("11GeV_R>0.4");
  TH2D * g0 = (TH2D *) fs->Get("8.8GeV_All");
  TH2D * g1 = (TH2D *) fs->Get("8.8GeV_R<0.4");
  TH2D * g2 = (TH2D *) fs->Get("8.8GeV_R>0.4");

  h0->Add(g0);
  h1->Add(g1);
  h2->Add(g2);

  double m = h0->GetMaximum();
  h1->SetMaximum(m);
  h2->SetMaximum(m);

  SetStyle(h0);
  SetStyle(h1);
  SetStyle(h2);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetLeftMargin(0.15);
  c0->SetBottomMargin(0.15);
  c0->SetLogz();

  h0->DrawClone("colz");

  c0->Print("zpt1.pdf");

  h1->DrawClone("colz");
  
  c0->Print("zpt2.pdf");

  h2->DrawClone("colz");

  c0->Print("zpt3.pdf");


  return 0;
}
  
  
