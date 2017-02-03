#include <iostream>
#include <fstream>
#include <cmath>

#include "TStyle.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLegend.h"

using namespace std;

int main(int argc, char * argv){
  gStyle->SetOptStat(0);

  double x, Q2, z, W, Estat;

  TChain * TCLAS12 = new TChain("data", "data");
  TCLAS12->Add("clas12_proton_pip.root");
  TCLAS12->Add("clas12_proton_pim.root");
  TCLAS12->SetBranchAddress("x", &x);
  TCLAS12->SetBranchAddress("Q2", &Q2);
  TCLAS12->SetBranchAddress("z", &z);
  TCLAS12->SetBranchAddress("Estat", &Estat);

  TChain * TSBS = new TChain("data", "data");
  TSBS->Add("sbs_neutron_pip_11.root");
  TSBS->Add("sbs_neutron_pim_11.root");
  TSBS->Add("sbs_neutron_pip_8.root");
  TSBS->Add("sbs_neutron_pim_8.root");
  TSBS->SetBranchAddress("x", &x);
  TSBS->SetBranchAddress("Q2", &Q2);
  TSBS->SetBranchAddress("z", &z);
  TSBS->SetBranchAddress("Estat", &Estat);

  TChain * TSOLIDp = new TChain("data", "data");
  TSOLIDp->Add("binsP11p.root");
  TSOLIDp->Add("binsP11m.root");
  TSOLIDp->Add("binsP8p.root");
  TSOLIDp->Add("binsP8m.root");
  TSOLIDp->SetBranchAddress("x", &x);
  TSOLIDp->SetBranchAddress("Q2", &Q2);
  TSOLIDp->SetBranchAddress("z", &z);
  TSOLIDp->SetBranchAddress("E0stat", &Estat);

  TChain * TSOLIDn = new TChain("data", "data");
  TSOLIDn->Add("binsN11p.root");
  TSOLIDn->Add("binsN11m.root");
  TSOLIDn->Add("binsN8p.root");
  TSOLIDn->Add("binsN8m.root");
  TSOLIDn->SetBranchAddress("x", &x);
  TSOLIDn->SetBranchAddress("Q2", &Q2);
  TSOLIDn->SetBranchAddress("z", &z);
  TSOLIDn->SetBranchAddress("E0stat", &Estat);


  double xSOLID[9] = {0.0, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 1.0};
  double xCLAS12[9] = {0.0, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.68, 1.0};
  double xSBS[9] = {0.0, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.78};

  double Q2SOLID[9] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};
  double Q2CLAS12[9] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.999, 6.0, 8.0, 10.0};
  double Q2SBS[9] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

  TH1D * hxSOLIDp = new TH1D("hxSOLIDp", "", 8, xSOLID);
  TH1D * hxSOLIDn = new TH1D("hxSOLIDn", "", 8, xSOLID);
  TH1D * hxCLAS12 = new TH1D("hxCLAS12", "", 8, xCLAS12);
  TH1D * hxSBS = new TH1D("hxSBS", "", 8, xSBS);
  
  TH1D * hQ2SOLIDp = new TH1D("hQ2SOLIDp", "", 8, Q2SOLID);
  TH1D * hQ2SOLIDn = new TH1D("hQ2SOLIDn", "", 8, Q2SOLID);
  TH1D * hQ2CLAS12 = new TH1D("hQ2CLAS12", "", 8, Q2CLAS12);
  TH1D * hQ2SBS = new TH1D("hQ2SBS", "", 8, Q2SBS);

  double Nt = 0;

  //CLAS12
  Nt = TCLAS12->GetEntries();
  for (int i = 0; i < Nt; i++){
    TCLAS12->GetEntry(i);
    if (z < 0.3 || z > 0.7) continue;
    W = sqrt(0.93827 * 0.93827 + Q2 * (1.0 / x - 1.0));
    //if (W < 2.3) continue;
    hxCLAS12->Fill(x, 1.0/Estat/Estat);
    hQ2CLAS12->Fill(Q2, 1.0/Estat/Estat);
  }
  for (int i = 1; i < 9; i++){
    hxCLAS12->SetBinContent(i, hxCLAS12->GetBinContent(i) / (xCLAS12[i] - xCLAS12[i-1]));
    hQ2CLAS12->SetBinContent(i, hQ2CLAS12->GetBinContent(i) / (Q2CLAS12[i] - Q2CLAS12[i-1]));
  }
  //hQ2CLAS12->SetBinContent(5, (hQ2CLAS12->GetBinContent(5) + hQ2CLAS12->GetBinContent(5)) / 2.0 * 2.0);
  //hQ2CLAS12->SetBinContent(6, hQ2CLAS12->GetBinContent(5) / 3.0);
  hxCLAS12->SetMarkerStyle(0);
  hxCLAS12->SetLineWidth(2);
  hxCLAS12->SetLineColor(2);
  hQ2CLAS12->SetMarkerStyle(0);
  hQ2CLAS12->SetLineWidth(2);
  hQ2CLAS12->SetLineColor(2);

  //SOLID p
  Nt = TSOLIDp->GetEntries();
  for (int i = 0; i < Nt; i++){
    TSOLIDp->GetEntry(i);
    if (z < 0.3 || z > 0.7) continue;
    W = sqrt(0.93827 * 0.93827 + Q2 * (1.0 / x - 1.0));
    if (W < 2.3) continue;
    Estat = Estat * 1.3;
    hxSOLIDp->Fill(x, 1.0/Estat/Estat);
    hQ2SOLIDp->Fill(Q2, 1.0/Estat/Estat);
  }
  for (int i = 1; i < 9; i++){
    hxSOLIDp->SetBinContent(i, hxSOLIDp->GetBinContent(i) / (xSOLID[i] - xSOLID[i-1]));
    hQ2SOLIDp->SetBinContent(i, hQ2SOLIDp->GetBinContent(i) / (Q2SOLID[i] - Q2SOLID[i-1]));
  }
  hxSOLIDp->SetMarkerStyle(0);
  hxSOLIDp->SetLineWidth(2);
  hxSOLIDp->SetLineColor(4);
  hQ2SOLIDp->SetMarkerStyle(0);
  hQ2SOLIDp->SetLineWidth(2);
  hQ2SOLIDp->SetLineColor(4);

  //SBS
  Nt = TSBS->GetEntries();
  for (int i = 0; i < Nt; i++){
    TSBS->GetEntry(i);
    if (z < 0.3 || z > 0.7) continue;
    W = sqrt(0.93827 * 0.93827 + Q2 * (1.0 / x - 1.0));
    if (W < 2.3) continue;
    hxSBS->Fill(x, 1.0/Estat/Estat);
    hQ2SBS->Fill(Q2, 1.0/Estat/Estat);
  }
  for (int i = 1; i < 9; i++){
    hxSBS->SetBinContent(i, hxSBS->GetBinContent(i) / (xSBS[i] - xSBS[i-1]));
    hQ2SBS->SetBinContent(i, hQ2SBS->GetBinContent(i) / (Q2SBS[i] - Q2SBS[i-1]));
  }
  hQ2SBS->SetBinContent(1, hQ2SBS->GetBinContent(2));
  hxSBS->SetBinContent(8, hxSBS->GetBinContent(7));
  hxSBS->SetMarkerStyle(0);
  hxSBS->SetLineWidth(2);
  hxSBS->SetLineColor(2);
  hQ2SBS->SetMarkerStyle(0);
  hQ2SBS->SetLineWidth(2);
  hQ2SBS->SetLineColor(2);

  //SOLID n
  Nt = TSOLIDn->GetEntries();
  for (int i = 0; i < Nt; i++){
    TSOLIDn->GetEntry(i);
    if (z < 0.3 || z > 0.7) continue;
    W = sqrt(0.93827 * 0.93827 + Q2 * (1.0 / x - 1.0));
    if (W < 2.3) continue;
    hxSOLIDn->Fill(x, 1.0/Estat/Estat);
    hQ2SOLIDn->Fill(Q2, 1.0/Estat/Estat);
  }
  for (int i = 0; i < 9; i++){
    hxSOLIDn->SetBinContent(i, hxSOLIDn->GetBinContent(i) / (xSOLID[i] - xSOLID[i-1]));
    hQ2SOLIDn->SetBinContent(i, hQ2SOLIDn->GetBinContent(i) / (Q2SOLID[i] - Q2SOLID[i-1]));
  }
  hxSOLIDn->SetMarkerStyle(0);
  hxSOLIDn->SetLineWidth(2);
  hxSOLIDn->SetLineColor(4);
  hQ2SOLIDn->SetMarkerStyle(0);
  hQ2SOLIDn->SetLineWidth(2);
  hQ2SOLIDn->SetLineColor(4);

  //Base
  TH1D * hx = new TH1D("hx", "", 1, 0.0, 0.8);
  hx->GetXaxis()->SetTitle("x");
  hx->GetXaxis()->CenterTitle(true);
  hx->GetXaxis()->SetTitleSize(0.06);
  hx->GetXaxis()->SetTitleOffset(1.15);
  hx->GetXaxis()->SetLabelSize(0.055);
  hx->GetYaxis()->SetTitle("FOM: (#deltaA_{UT})^{-2} / #Deltax");
  hx->GetYaxis()->CenterTitle(true);
  hx->GetYaxis()->SetTitleSize(0.06);
  hx->GetYaxis()->SetTitleOffset(1.15);
  hx->GetYaxis()->SetLabelSize(0.055);
  hx->GetYaxis()->SetRangeUser(1.e3, 5e9);

  TH1D * hQ2 = new TH1D("hQ2", "", 1, 0.0, 11.0);
  hQ2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  hQ2->GetXaxis()->CenterTitle(true);
  hQ2->GetXaxis()->SetTitleSize(0.06);
  hQ2->GetXaxis()->SetTitleOffset(1.15);
  hQ2->GetXaxis()->SetLabelSize(0.055);
  hQ2->GetYaxis()->SetTitle("FOM: (#deltaA_{UT})^{-2} / #DeltaQ^{2} (GeV^{-2})");
  hQ2->GetYaxis()->CenterTitle(true);
  hQ2->GetYaxis()->SetTitleSize(0.06);
  hQ2->GetYaxis()->SetTitleOffset(1.15);
  hQ2->GetYaxis()->SetLabelSize(0.055);
  hQ2->GetYaxis()->SetRangeUser(1.e2, 5e8);

  
  TLegend * leg_sbs = new TLegend(0.15, 0.15, 0.75, 0.28);
  leg_sbs->AddEntry(hxSOLIDn, "SoLID E12-10-006 with ^{3}He", "l");
  leg_sbs->AddEntry(hxSBS, "SBS E12-09-018 with ^{3}He", "l");

  TLegend * leg_clas = new TLegend(0.15, 0.15, 0.75, 0.28);
  leg_clas->AddEntry(hxSOLIDp, "SoLID E12-11-108 with NH_{3}", "l");
  leg_clas->AddEntry(hxCLAS12, "CLAS12 C12-11-111 with HDice", "l");


  TCanvas * cx = new TCanvas("cx", "", 800, 600);
  cx->SetBottomMargin(0.15);
  cx->SetLeftMargin(0.15);
  cx->SetLogy();
  hx->Draw();
  hxSOLIDn->Draw("esame");
  hxSBS->Draw("esame");
  leg_sbs->Draw("same");
  cx->Print("fom_solid_sbs_x.pdf");
  hx->Draw();
  hxSOLIDp->Draw("esame");
  hxCLAS12->Draw("esame");
  leg_clas->Draw("same");
  cx->Print("fom_solid_clas12_x.pdf");

  TCanvas * cQ2 = new TCanvas("cQ2", "", 800, 600);
  cQ2->SetBottomMargin(0.15);
  cQ2->SetLeftMargin(0.15);
  cQ2->SetLogy();
  hQ2->Draw();
  hQ2SOLIDn->Draw("esame");
  hQ2SBS->Draw("esame");
  leg_sbs->Draw("same");
  cQ2->Print("fom_solid_sbs_Q2.pdf");
  hQ2->Draw();
  hQ2SOLIDp->Draw("esame");
  hQ2CLAS12->Draw("esame");
  leg_clas->Draw("same");
  cQ2->Print("fom_solid_clas12_Q2.pdf");


  return 0;
}
  

  
				 

