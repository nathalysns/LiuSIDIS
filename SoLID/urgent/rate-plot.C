#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./rate-plot <rootfile>" << endl;
    return 0;
  }

  gStyle->SetOptStat(0);
  
  TFile * fs = new TFile(argv[1], "r");
  TH2D * hxQ2 = (TH2D *) fs->Get("hxQ2");
  TH2D * h0xQ2 = (TH2D *) fs->Get("h0xQ2");
  TH2D * h1xQ2 = (TH2D *) fs->Get("h1xQ2");
  TH2D * h2xQ2 = (TH2D *) fs->Get("h2xQ2");
  TH2D * h3xQ2 = (TH2D *) fs->Get("h3xQ2");
  TH2D * hzPt = (TH2D *) fs->Get("hzPt");
  TH2D * h0zPt = (TH2D *) fs->Get("h0zPt");
  TH2D * h1zPt = (TH2D *) fs->Get("h1zPt");
  TH2D * h2zPt = (TH2D *) fs->Get("h2zPt");
  TH2D * h3zPt = (TH2D *) fs->Get("h3zPt");

  h0xQ2->Divide(hxQ2);
  h1xQ2->Divide(hxQ2);
  h2xQ2->Divide(hxQ2);
  h3xQ2->Divide(hxQ2);
  h0zPt->Divide(hzPt);
  h1zPt->Divide(hzPt);
  h2zPt->Divide(hzPt);
  h3zPt->Divide(hzPt);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetBottomMargin(0.15);
  c0->SetLeftMargin(0.15);

  TH2D * hBxQ2 = new TH2D("hBxQ2", "", 1, 0.0, 1.0, 1, 0.0, 10.0);
  hBxQ2->GetXaxis()->SetTitle("x");
  hBxQ2->GetXaxis()->SetLabelSize(0.055);
  hBxQ2->GetXaxis()->CenterTitle(true);
  hBxQ2->GetXaxis()->SetTitleSize(0.055);
  hBxQ2->GetXaxis()->SetTitleOffset(1.15);
  hBxQ2->GetYaxis()->SetTitle("Q^{2}(GeV^{2})");
  hBxQ2->GetYaxis()->SetLabelSize(0.055);
  hBxQ2->GetYaxis()->CenterTitle(true);
  hBxQ2->GetYaxis()->SetTitleSize(0.055);
  hBxQ2->GetYaxis()->SetTitleOffset(1.15);
  hBxQ2->SetMinimum(0);
  hBxQ2->SetMaximum(1.01);


  hBxQ2->DrawClone("axis");
  h0xQ2->DrawClone("colzsame");
  c0->Print("rateplotxQ2.pdf(");

  hBxQ2->DrawClone("axis");
  h1xQ2->DrawClone("colzsame");
  c0->Print("rateplotxQ2.pdf");

  hBxQ2->DrawClone("axis");
  h2xQ2->DrawClone("colzsame");
  c0->Print("rateplotxQ2.pdf");

  hBxQ2->DrawClone("axis");
  h3xQ2->DrawClone("colzsame");
  c0->Print("rateplotxQ2.pdf)");

  TH2D * hBzPt = new TH2D("hBzPt", "", 1, 0.0, 1.0, 1, 0.0, 2.0);
  hBzPt->GetXaxis()->SetTitle("z");
  hBzPt->GetXaxis()->SetLabelSize(0.055);
  hBzPt->GetXaxis()->CenterTitle(true);
  hBzPt->GetXaxis()->SetTitleSize(0.055);
  hBzPt->GetXaxis()->SetTitleOffset(1.15);
  hBzPt->GetYaxis()->SetTitle("P_{hT}(GeV)");
  hBzPt->GetYaxis()->SetLabelSize(0.055);
  hBzPt->GetYaxis()->CenterTitle(true);
  hBzPt->GetYaxis()->SetTitleSize(0.055);
  hBzPt->GetYaxis()->SetTitleOffset(1.15);
  hBzPt->SetMinimum(0);
  hBzPt->SetMaximum(1);

  hBzPt->DrawClone("axis");
  h0zPt->DrawClone("colzsame");
  c0->Print("rateplotzPt.pdf(");

  hBzPt->DrawClone("axis");
  h1zPt->DrawClone("colzsame");
  c0->Print("rateplotzPt.pdf");

  hBzPt->DrawClone("axis");
  h2zPt->DrawClone("colzsame");
  c0->Print("rateplotzPt.pdf");
  
  hBzPt->DrawClone("axis");
  h3zPt->DrawClone("colzsame");
  c0->Print("rateplotzPt.pdf)");
  

  return 0;
}

  
