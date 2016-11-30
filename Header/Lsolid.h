#ifndef _LSOLID_H_
#define _LSOLID_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

class Lsolid_3he{
 protected:
  TFile * file_negative;
  TFile * file_positive;
  TH2F * acc_FA_e;
  TH2F * acc_LA_e;
  TH2F * acc_FA_pip;
  TH2F * acc_FA_pim;
 public:
  Lsolid_3he();
  int Initialize();//get acceptance files
  double GetAcceptance(const TLorentzVector p, const char * option);//get acceptance
};

Lsolid_3he::Lsolid_3he(){
}

int Lsolid_3he::Initialize(){
  file_negative = new TFile("../Header/SoLIDfiles/acceptance_solid_CLEO_SIDIS_3he_negative_output.root","r");
  file_positive = new TFile("../Header/SoLIDfiles/acceptance_solid_CLEO_SIDIS_3he_positive_output.root","r");
  acc_FA_e = (TH2F *) file_negative->Get("acceptance_forwardangle");
  acc_LA_e = (TH2F *) file_negative->Get("acceptance_largeangle");
  acc_FA_pip = (TH2F *) file_positive->Get("acceptance_forwardangle");
  acc_FA_pim = (TH2F *) file_negative->Get("acceptance_forwardangle");
  return 0;
}

double Lsolid_3he::GetAcceptance(const TLorentzVector p, const char * option){
  double theta = p.Theta() / M_PI * 180.0;
  double mom = p.P();
  if (theta < 8.0 || theta > 30.0) return 0;
  double acc = 0.0;
  if (strcmp(option, "e-") == 0){
    acc = acc_FA_e->GetBinContent(acc_FA_e->GetXaxis()->FindBin(theta), acc_FA_e->GetYaxis()->FindBin(mom));
    if (mom > 3.5){
      acc += acc_LA_e->GetBinContent(acc_LA_e->GetXaxis()->FindBin(theta), acc_LA_e->GetYaxis()->FindBin(mom));
    }
  }
  else if (strcmp(option, "pi+") == 0){
    acc = acc_FA_pip->GetBinContent(acc_FA_pip->GetXaxis()->FindBin(theta), acc_FA_pip->GetYaxis()->FindBin(mom));
  }
  else if (strcmp(option, "pi-") == 0){
    acc = acc_FA_pim->GetBinContent(acc_FA_pim->GetXaxis()->FindBin(theta), acc_FA_pim->GetYaxis()->FindBin(mom));
  }
  return acc;
}









#endif
