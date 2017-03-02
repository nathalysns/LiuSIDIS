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
  
  TString file = argv[1];
  TString file1 = file + ".root";
  TString file2 = file + "_R0.4.root";
  string obj = argv[2];
  TString save = argv[3];
  TString save1 = save + ".png";
  TString save2 = save + "_R0.4.png";
 
  TFile * fA = new TFile(file1.Data(), "r");
  TH2D * FA = (TH2D *) fA->Get(obj.data());
  TFile * fR = new TFile(file2.Data(), "r");
  TH2D * FR = (TH2D *) fR->Get(obj.data());

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

  FR->GetXaxis()->CenterTitle(true);
  FR->GetXaxis()->SetTitleFont(42);
  FR->GetXaxis()->SetTitleSize(0.06);
  FR->GetXaxis()->SetTitleOffset(1.15);
  FR->GetXaxis()->SetLabelFont(42);
  FR->GetXaxis()->SetLabelSize(0.055);
  FR->GetYaxis()->CenterTitle(true);
  FR->GetYaxis()->SetTitleFont(42);
  FR->GetYaxis()->SetTitleSize(0.06);
  FR->GetYaxis()->SetTitleOffset(1.15);
  FR->GetYaxis()->SetLabelFont(42);
  FR->GetYaxis()->SetLabelSize(0.055);

  double zmax = FA->GetMaximum() * 1.01;
  double zmin = zmax / 10000.;
  FA->GetZaxis()->SetRangeUser(zmin, zmax);
  FR->GetZaxis()->SetRangeUser(zmin, zmax);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetLogz();
  c0->SetLeftMargin(0.15);
  c0->SetBottomMargin(0.15);
  FA->Draw("colz");
  c0->Print(save1.Data());

  TCanvas * c1 = new TCanvas("c1", "", 800, 600);
  c1->SetLogz();
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  FR->Draw("colz");
  c1->Print(save2.Data());

  return 0;
}
