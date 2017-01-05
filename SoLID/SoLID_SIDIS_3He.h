#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"

#include "Lsidis.h"

// Acceptance 
TFile * file_e = new TFile("Acceptance/acceptance_solid_SIDIS_He3_electron_201701_1e7_output.root", "r");
TFile * file_pi = new TFile("Acceptance/acceptance_solid_SIDIS_He3_pim_201701_1e7_output.root", "r");
TH2F * acc_FA_e = (TH2F *) file_e->Get("acceptance_ThetaP_forwardangle");
TH2F * acc_LA_e = (TH2F *) file_e->Get("acceptance_ThetaP_largeangle");
TH2F * acc_FA_pi = (TH2F *) file_pi->Get("acceptance_ThetaP_forwardangle");
TH2F * acc_LA_pi = (TH2F *) file_pi->Get("acceptance_ThetaP_largeangle");

double GetAcceptance_e(const TLorentzVector p){//Get electron acceptance
  double theta = p.Theta() / M_PI * 180.0;
  if (theta < 8.0 || theta > 30.0) return 0;
  double mom = p.P();
  double acc = 0;
  acc += acc_FA_e->GetBinContent(acc_FA_e->GetXaxis()->FindBin(theta), acc_FA_e->GetYaxis()->FindBin(mom));
  if (mom > 3.5)
    acc += acc_LA_e->GetBinContent(acc_LA_e->GetXaxis()->FindBin(theta), acc_LA_e->GetYaxis()->FindBin(mom));
  return acc;
}

double GetAcceptance_pi(const TLorentzVector p){//Get pion acceptance
  double theta = p.Theta() / M_PI * 180.0;
  if (theta < 8.0 || theta > 18.0) return 0;
  double mom = p.P();
  double acc = 0;
  acc += acc_FA_pi->GetBinContent(acc_FA_pi->GetXaxis()->FindBin(theta), acc_FA_pi->GetYaxis()->FindBin(mom));
  return acc;
}
