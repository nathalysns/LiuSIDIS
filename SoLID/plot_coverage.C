#include <iostream>
#include <fstream>
#include <cstring>

#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;

int main(int argc, char * argv[]){

  if (argc < 4){
    cout << "./plot <rootfile> <object> <saveplot>" << endl;
    return 1;
  }

  string obj = argv[2];
  string obj_FA = obj + "_FA";
  string obj_LA = obj + "_LA";

  TFile * fs = new TFile(argv[1], "r");
  TH2D * FA = (TH2D *) fs->Get(obj_FA.data());
  TH2D * LA = (TH2D *) fs->Get(obj_LA.data());

  FA->SetMarkerColor(1);
  FA->GetXaxis()->CenterTitle(true);
  FA->GetXaxis()->SetTitleFont(42);
  FA->GetXaxis()->SetTitleSize(0.06);
  FA->GetXaxis()->SetTitleOffset(1.15);
  FA->GetXaxis()->SetLabelFont(42);
  FA->GetXaxis()->SetLabelSize(0.055);
  FA->GetYaxis()->CenterTitle(true);
  FA->GetYaxis()->SetTitleFont(42);
  FA->GetYaxis()->SetTitleSize(0.06);
  FA->GetYaxis()->SetTitleOffset(1.15);
  FA->GetYaxis()->SetLabelFont(42);
  FA->GetYaxis()->SetLabelSize(0.055);

  LA->SetMarkerColor(3);
  LA->GetXaxis()->CenterTitle(true);
  LA->GetXaxis()->SetTitleFont(42);
  LA->GetXaxis()->SetTitleSize(0.06);
  LA->GetXaxis()->SetTitleOffset(1.15);
  LA->GetXaxis()->SetLabelFont(42);
  LA->GetXaxis()->SetLabelSize(0.055);
  LA->GetYaxis()->CenterTitle(true);
  LA->GetYaxis()->SetTitleFont(42);
  LA->GetYaxis()->SetTitleSize(0.06);
  LA->GetYaxis()->SetTitleOffset(1.15);
  LA->GetYaxis()->SetLabelFont(42);
  LA->GetYaxis()->SetLabelSize(0.055);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetLeftMargin(0.15);
  c0->SetBottomMargin(0.15);
  FA->Draw();
  LA->Draw("same");

  c0->Print(argv[3]);

  return 0;
}
