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
 
  TFile * fs = new TFile(argv[1], "r");
  TH2D * FA = (TH2D *) fs->Get(obj.data());

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

  double zmax = FA->GetMaximum() * 1.01;
  double zmin = zmax / 10000.;
  FA->GetZaxis()->SetRangeUser(zmin, zmax);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetLogz();
  c0->SetLeftMargin(0.15);
  c0->SetBottomMargin(0.15);
  FA->Draw("colz");

  c0->Print(argv[3]);

  return 0;
}
