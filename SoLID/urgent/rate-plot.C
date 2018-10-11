#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./rate-plot <rootfile> <prename>" << endl;
    return 0;
  }

  TString pre = argv[2];
  gStyle->SetOptStat(0);
  
  TFile * fs = new TFile(argv[1], "r");
  TH2D * h1xQ2 = (TH2D *) fs->Get("h1xQ2");
  TH2D * h2xQ2 = (TH2D *) fs->Get("h2xQ2");
  TH2D * h1zPt = (TH2D *) fs->Get("h1zPt");
  TH2D * h2zPt = (TH2D *) fs->Get("h2zPt");

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
  h2xQ2->Divide(h1xQ2);
  h2xQ2->DrawClone("colzsame");

  TString filexQ2 = pre + "xQ2.pdf";
  c0->Print(filexQ2.Data());

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
  h2zPt->Divide(h1zPt);
  h2zPt->DrawClone("colzsame");

  TString filezPt = pre + "zPt.pdf";
  c0->Print(filezPt);

  return 0;
}

  
