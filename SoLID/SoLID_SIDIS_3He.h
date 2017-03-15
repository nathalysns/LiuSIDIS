#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TChain.h"

#include "Lsidis.h"

// Acceptance 
TFile * file_e = new TFile("Acceptance/acceptance_solid_SIDIS_He3_electron_201701_1e7_output.root", "r");
TFile * file_pi = new TFile("Acceptance/acceptance_solid_SIDIS_He3_pim_201701_1e7_output.root", "r");
//TFile * file_e = new TFile("Acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root", "r");
//TFile * file_pi = new TFile("Acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root", "r");
TH2F * acc_FA_e = (TH2F *) file_e->Get("acceptance_ThetaP_forwardangle");
TH2F * acc_LA_e = (TH2F *) file_e->Get("acceptance_ThetaP_largeangle");
TH2F * acc_FA_pi = (TH2F *) file_pi->Get("acceptance_ThetaP_forwardangle");
TH2F * acc_LA_pi = (TH2F *) file_pi->Get("acceptance_ThetaP_largeangle");
//TH2F * acc_FA_e = (TH2F *) file_e->Get("acceptance_forwardangle");
//TH2F * acc_LA_e = (TH2F *) file_e->Get("acceptance_largeangle");
//TH2F * acc_FA_pi = (TH2F *) file_pi->Get("acceptance_forwardangle");
//TH2F * acc_LA_pi = (TH2F *) file_pi->Get("acceptance_largeangle");


double Rfactor0 = 1.0e5;

double GetAcceptance_e(const TLorentzVector p, const char * detector = "all"){//Get electron acceptance
  double theta = p.Theta() / M_PI * 180.0;
  if (theta < 8.0 || theta > 30.0) return 0;
  double mom = p.P();
  double acc = 0;
  if (strcmp(detector, "FA") == 0 || strcmp(detector, "all") == 0)
    acc += acc_FA_e->GetBinContent(acc_FA_e->GetXaxis()->FindBin(theta), acc_FA_e->GetYaxis()->FindBin(mom));
  if (mom > 3.5 && (strcmp(detector, "LA") == 0 || strcmp(detector, "all") == 0))
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

int GetTotalRate(const double Ebeam, const char * hadron){//Estimate the total rate
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2, 1);
  sidis.SetHadron(hadron);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+10 * pow(0.197327, 2);
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 10.0, 0.7, 1.8, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  double sum = 0.0;
  Long64_t Nsim = 100000000;
  double weight = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%(Nsim/20) == 0) std::cout << i << std::endl;
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      sidis.CalculateRfactor();
      if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      sum += weight * GetAcceptance_e(lp) * GetAcceptance_pi(Ph);
    }
  }
  printf("\n");
  printf("Total rate: %.4E  (%.1f GeV %s)\n\n", sum * lumi / Nsim, Ebeam, hadron);
  return 0;
}

int MakeKinematicCoveragePlots(const double Ebeam, const char * savefile){
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2,1);
  sidis.SetHadron("pi+");
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 9.0, 0.7, 2.0, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  TFile * fs = new TFile(savefile, "RECREATE");
  gStyle->SetOptStat(0);
  //(x, Q2)
  TH2D * xQ2_FA = new TH2D("xQ2_FA", "", 700, 0.0, 0.7, 900, 0.0, 9.0);
  xQ2_FA->GetXaxis()->SetTitle("x");
  xQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  TH2D * xQ2_LA = new TH2D("xQ2_LA", "", 700, 0.0, 0.7, 900, 0.0, 9.0);
  xQ2_LA->GetXaxis()->SetTitle("x");
  xQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(x, W)
  TH2D * xW_FA = new TH2D("xW_FA", "", 700, 0.0, 0.7, 500, 2.0, 4.5);
  xW_FA->GetXaxis()->SetTitle("x");
  xW_FA->GetYaxis()->SetTitle("W / GeV");
  TH2D * xW_LA = new TH2D("xW_LA", "", 700, 0.0, 0.7, 500, 2.0, 4.5);
  xW_LA->GetXaxis()->SetTitle("x");
  xW_LA->GetYaxis()->SetTitle("W / GeV");
  //(x, Wp)
  TH2D * xWp_FA = new TH2D("xWp_FA", "", 700, 0.0, 0.7, 500, 1.5, 4.0);
  xWp_FA->GetXaxis()->SetTitle("x");
  xWp_FA->GetYaxis()->SetTitle("W' / GeV");
  TH2D * xWp_LA = new TH2D("xWp_LA", "", 700, 0.0, 0.7, 500, 1.5, 4.0);
  xWp_LA->GetXaxis()->SetTitle("x");
  xWp_LA->GetYaxis()->SetTitle("W' / GeV");
  //(x, z)
  TH2D * xz_FA = new TH2D("xz_FA", "", 700, 0.0, 0.7, 600, 0.2, 0.8);
  xz_FA->GetXaxis()->SetTitle("x");
  xz_FA->GetYaxis()->SetTitle("z");
  TH2D * xz_LA = new TH2D("xz_LA", "", 700, 0.0, 0.7, 600, 0.2, 0.8);
  xz_LA->GetXaxis()->SetTitle("x");
  xz_LA->GetYaxis()->SetTitle("z");
  //(x, Pt)
  TH2D * xPt_FA = new TH2D("xPt_FA", "", 700, 0.0, 0.7, 800, 0.0, 2.0);
  xPt_FA->GetXaxis()->SetTitle("x");
  xPt_FA->GetYaxis()->SetTitle("P_{T} / GeV");
  TH2D * xPt_LA = new TH2D("xPt_LA", "", 700, 0.0, 0.7, 800, 0.0, 2.0);
  xPt_LA->GetXaxis()->SetTitle("x");
  xPt_LA->GetYaxis()->SetTitle("P_{T} / GeV");
  //(z, Pt)
  TH2D * zPt_FA = new TH2D("zPt_FA", "", 600, 0.2, 0.8, 800, 0.0, 2.0);
  zPt_FA->GetXaxis()->SetTitle("z");
  zPt_FA->GetYaxis()->SetTitle("P_{T} / GeV");
  TH2D * zPt_LA = new TH2D("zPt_LA", "", 600, 0.2, 0.8, 800, 0.0, 2.0);
  zPt_LA->GetXaxis()->SetTitle("z");
  zPt_LA->GetYaxis()->SetTitle("P_{T} / GeV");
  //(z, Q2)
  TH2D * zQ2_FA = new TH2D("zQ2_FA", "", 600, 0.2, 0.8, 900, 0.0, 9.0);
  zQ2_FA->GetXaxis()->SetTitle("z");
  zQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  TH2D * zQ2_LA = new TH2D("zQ2_LA", "", 600, 0.2, 0.8, 900, 0.0, 9.0);
  zQ2_LA->GetXaxis()->SetTitle("z");
  zQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(z, W)
  TH2D * zW_FA = new TH2D("zW_FA", "", 600, 0.2, 0.8, 500, 2.0, 4.5);
  zW_FA->GetXaxis()->SetTitle("z");
  zW_FA->GetYaxis()->SetTitle("W / GeV");
  TH2D * zW_LA = new TH2D("zW_LA", "", 600, 0.2, 0.8, 500, 2.0, 4.5);
  zW_LA->GetXaxis()->SetTitle("z");
  zW_LA->GetYaxis()->SetTitle("W / GeV");
  //(z, Wp)
  TH2D * zWp_FA = new TH2D("zWp_FA", "", 600, 0.2, 0.8, 500, 1.5, 4.0);
  zWp_FA->GetXaxis()->SetTitle("z");
  zWp_FA->GetYaxis()->SetTitle("W' / GeV");
  TH2D * zWp_LA = new TH2D("zWp_LA", "", 600, 0.2, 0.8, 500, 1.5, 4.0);
  zWp_LA->GetXaxis()->SetTitle("z");
  zWp_LA->GetYaxis()->SetTitle("W' / GeV");
  //(Pt, Q2)
  TH2D * PtQ2_FA = new TH2D("PtQ2_FA", "", 800, 0.0, 2.0, 900, 0.0, 9.0);
  PtQ2_FA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  TH2D * PtQ2_LA = new TH2D("PtQ2_LA", "", 800, 0.0, 2.0, 900, 0.0, 9.0);
  PtQ2_LA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(Pt, W)
  TH2D * PtW_FA = new TH2D("PtW_FA", "", 800, 0.0, 2.0, 500, 2.0, 4.5);
  PtW_FA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtW_FA->GetYaxis()->SetTitle("W / GeV");
  TH2D * PtW_LA = new TH2D("PtW_LA", "", 800, 0.0, 2.0, 500, 2.0, 4.5);
  PtW_LA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtW_LA->GetYaxis()->SetTitle("W / GeV");
  //(Pt, Wp)
  TH2D * PtWp_FA = new TH2D("PtWp_FA", "", 800, 0.0, 2.0, 500, 1.5, 4.0);
  PtWp_FA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtWp_FA->GetYaxis()->SetTitle("W' / GeV");
  TH2D * PtWp_LA = new TH2D("PtWp_LA", "", 800, 0.0, 2.0, 500, 1.5, 4.0);
  PtWp_LA->GetXaxis()->SetTitle("P_{T} / GeV");
  PtWp_LA->GetYaxis()->SetTitle("W' / GeV");
  //(W, Q2)
  TH2D * WQ2_FA = new TH2D("WQ2_FA", "", 500, 2.0, 4.5, 900, 0.0, 9.0);
  WQ2_FA->GetXaxis()->SetTitle("W / GeV");
  WQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  TH2D * WQ2_LA = new TH2D("WQ2_LA", "", 500, 2.0, 4.5, 900, 0.0, 9.0);
  WQ2_LA->GetXaxis()->SetTitle("W / GeV");
  WQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(Wp, Q2)
  TH2D * WpQ2_FA = new TH2D("WpQ2_FA", "", 500, 1.5, 4.0, 900, 0.0, 9.0);
  WpQ2_FA->GetXaxis()->SetTitle("W' / GeV");
  WpQ2_FA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  TH2D * WpQ2_LA = new TH2D("WpQ2_LA", "", 500, 1.5, 4.0, 900, 0.0, 9.0);
  WpQ2_LA->GetXaxis()->SetTitle("W' / GeV");
  WpQ2_LA->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  double x, Q2, z, Pt, W, Wp;
  double weight, acc_FA, acc_LA;
  TLorentzVector lp, Ph;
  for (Long64_t i = 0; i < 100000000; i++){
    if (i % 1000000 == 0) std::cout << i << " %" << std::endl;
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      z = sidis.GetVariable("z");
      if (z < 0.3 || z > 0.7) continue;
      Q2 = sidis.GetVariable("Q2");
      if (Q2 < 1.0) continue;
      W = sidis.GetVariable("W");
      if (W < 2.3) continue;
      Wp = sidis.GetVariable("Wp");
      if (Wp < 1.6) continue;
      x = sidis.GetVariable("x");
      Pt = sidis.GetVariable("Pt");
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc_FA = GetAcceptance_e(lp, "FA") * GetAcceptance_pi(Ph);
      acc_LA = GetAcceptance_e(lp, "LA") * GetAcceptance_pi(Ph);
      if (acc_FA > 0){
	xQ2_FA->Fill(x, Q2, acc_FA);
	xW_FA->Fill(x, W, acc_FA);
	xz_FA->Fill(x, z, acc_FA);
	xPt_FA->Fill(x, Pt, acc_FA);
	xWp_FA->Fill(x, Wp, acc_FA);
	zPt_FA->Fill(z, Pt, acc_FA);
	zQ2_FA->Fill(z, Q2, acc_FA);
	zW_FA->Fill(z, W, acc_FA);
	zWp_FA->Fill(z, Wp, acc_FA);
	PtQ2_FA->Fill(Pt, Q2, acc_FA);
	PtW_FA->Fill(Pt, W, acc_FA);
	PtWp_FA->Fill(Pt, Wp, acc_FA);
	WQ2_FA->Fill(W, Q2, acc_FA);
	WpQ2_FA->Fill(Wp, Q2, acc_FA);
      }
      if (acc_LA > 0){
	xQ2_LA->Fill(x, Q2, acc_LA);
	xW_LA->Fill(x, W, acc_LA);
	xz_LA->Fill(x, z, acc_LA);
	xPt_LA->Fill(x, Pt, acc_LA);
	xWp_LA->Fill(x, Wp, acc_LA);
	zPt_LA->Fill(z, Pt, acc_LA);
	zQ2_LA->Fill(z, Q2, acc_LA);
	zW_LA->Fill(z, W, acc_LA);
	zWp_LA->Fill(z, Wp, acc_LA);
	PtQ2_LA->Fill(Pt, Q2, acc_LA);
	PtW_LA->Fill(Pt, W, acc_LA);
	PtWp_LA->Fill(Pt, Wp, acc_LA);
	WQ2_LA->Fill(W, Q2, acc_LA);
	WpQ2_LA->Fill(Wp, Q2, acc_LA);
      }
    }
  }
  xQ2_FA->Divide(xQ2_FA); xQ2_FA->Scale(100);
  xQ2_LA->Divide(xQ2_LA); xQ2_LA->Scale(100);
  xW_FA->Divide(xW_FA); xW_FA->Scale(100);
  xW_LA->Divide(xW_LA); xW_LA->Scale(100);
  xz_FA->Divide(xz_FA); xz_FA->Scale(100);
  xz_LA->Divide(xz_LA); xz_LA->Scale(100);
  xPt_FA->Divide(xPt_FA); xPt_FA->Scale(100);
  xPt_LA->Divide(xPt_LA); xPt_LA->Scale(100);
  xWp_FA->Divide(xWp_FA); xWp_FA->Scale(100);
  xWp_LA->Divide(xWp_LA); xWp_LA->Scale(100);
  zW_FA->Divide(zW_FA); zW_FA->Scale(100);
  zW_LA->Divide(zW_LA); zW_LA->Scale(100);
  zQ2_FA->Divide(zQ2_FA); zQ2_FA->Scale(100);
  zQ2_LA->Divide(zQ2_LA); zQ2_LA->Scale(100);
  zWp_FA->Divide(zWp_FA); zWp_FA->Scale(100);
  zWp_LA->Divide(zWp_LA); zWp_LA->Scale(100);
  zPt_FA->Divide(zPt_FA); zPt_FA->Scale(100);
  zPt_LA->Divide(zPt_LA); zPt_LA->Scale(100);
  PtQ2_FA->Divide(PtQ2_FA); PtQ2_FA->Scale(100);
  PtQ2_LA->Divide(PtQ2_LA); PtQ2_LA->Scale(100);
  WQ2_FA->Divide(WQ2_FA); WQ2_FA->Scale(100);
  WQ2_LA->Divide(WQ2_LA); WQ2_LA->Scale(100);
  PtW_FA->Divide(PtW_FA); PtW_FA->Scale(100);
  PtW_LA->Divide(PtW_LA); PtW_LA->Scale(100);
  PtWp_FA->Divide(PtWp_FA); PtWp_FA->Scale(100);
  PtWp_LA->Divide(PtWp_LA); PtWp_LA->Scale(100);
  WpQ2_FA->Divide(WpQ2_FA); WpQ2_FA->Scale(100);
  WpQ2_LA->Divide(WpQ2_LA); WpQ2_LA->Scale(100);
  fs->Write();
  return 0;
}

int MakeRateDistributionPlots(const double Ebeam, const char * savefile){
  double lumi = 1.0e+10 * pow(0.197327, 2);
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2,1);
  sidis.SetHadron("pi+");
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 10.0, 0.7, 2.0, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  TFile * fs = new TFile(savefile, "RECREATE");
  gStyle->SetOptStat(0);
  //(x, Q2)
  TH2D * xQ2 = new TH2D("xQ2", "", 700, 0.0, 0.7, 900, 0.0, 9.0);
  xQ2->GetXaxis()->SetTitle("x");
  xQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(x, W)
  TH2D * xW = new TH2D("xW", "", 700, 0.0, 0.7, 500, 2.0, 4.5);
  xW->GetXaxis()->SetTitle("x");
  xW->GetYaxis()->SetTitle("W / GeV");
  //(x, Wp)
  TH2D * xWp = new TH2D("xWp", "", 700, 0.0, 0.7, 500, 1.5, 4.0);
  xWp->GetXaxis()->SetTitle("x");
  xWp->GetYaxis()->SetTitle("W' / GeV");
  //(x, z)
  TH2D * xz = new TH2D("xz", "", 700, 0.0, 0.7, 600, 0.2, 0.8);
  xz->GetXaxis()->SetTitle("x");
  xz->GetYaxis()->SetTitle("z");
  //(x, Pt)
  TH2D * xPt = new TH2D("xPt", "", 700, 0.0, 0.7, 800, 0.0, 2.0);
  xPt->GetXaxis()->SetTitle("x");
  xPt->GetYaxis()->SetTitle("P_{T} / GeV");
  //(z, Pt)
  TH2D * zPt = new TH2D("zPt", "", 600, 0.2, 0.8, 800, 0.0, 2.0);
  zPt->GetXaxis()->SetTitle("z");
  zPt->GetYaxis()->SetTitle("P_{T} / GeV");
  //(z, Q2)
  TH2D * zQ2 = new TH2D("zQ2", "", 600, 0.2, 0.8, 900, 0.0, 9.0);
  zQ2->GetXaxis()->SetTitle("z");
  zQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(z, W)
  TH2D * zW = new TH2D("zW", "", 600, 0.2, 0.8, 500, 2.0, 4.5);
  zW->GetXaxis()->SetTitle("z");
  zW->GetYaxis()->SetTitle("W / GeV");
  //(z, Wp)
  TH2D * zWp = new TH2D("zWp", "", 600, 0.2, 0.8, 500, 1.5, 4.0);
  zWp->GetXaxis()->SetTitle("z");
  zWp->GetYaxis()->SetTitle("W' / GeV");
  //(Pt, Q2)
  TH2D * PtQ2 = new TH2D("PtQ2", "", 800, 0.0, 2.0, 900, 0.0, 9.0);
  PtQ2->GetXaxis()->SetTitle("P_{T} / GeV");
  PtQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(Pt, W)
  TH2D * PtW = new TH2D("PtW", "", 800, 0.0, 2.0, 500, 2.0, 4.5);
  PtW->GetXaxis()->SetTitle("P_{T} / GeV");
  PtW->GetYaxis()->SetTitle("W / GeV");
  //(Pt, Wp)
  TH2D * PtWp = new TH2D("PtWp", "", 800, 0.0, 2.0, 500, 1.5, 4.0);
  PtWp->GetXaxis()->SetTitle("P_{T} / GeV");
  PtWp->GetYaxis()->SetTitle("W' / GeV");
  //(W, Q2)
  TH2D * WQ2 = new TH2D("WQ2", "", 500, 2.0, 4.5, 900, 0.0, 9.0);
  WQ2->GetXaxis()->SetTitle("W / GeV");
  WQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  //(Wp, Q2)
  TH2D * WpQ2 = new TH2D("WpQ2", "", 500, 1.5, 4.0, 900, 0.0, 9.0);
  WpQ2->GetXaxis()->SetTitle("W' / GeV");
  WpQ2->GetYaxis()->SetTitle("Q^{2} / GeV^{2}");
  double x, Q2, z, Pt, W, Wp;
  double weight, acc;
  TLorentzVector lp, Ph;
  Long64_t Nsim = 100000000;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i % 10000000 == 0) std::cout << i << " %" << std::endl;
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      z = sidis.GetVariable("z");
      if (z < 0.3 || z > 0.7) continue;
      Q2 = sidis.GetVariable("Q2");
      if (Q2 < 1.0) continue;
      W = sidis.GetVariable("W");
      if (W < 2.3) continue;
      Wp = sidis.GetVariable("Wp");
      if (Wp < 1.6) continue;
      x = sidis.GetVariable("x");
      Pt = sidis.GetVariable("Pt");
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp, "all") * GetAcceptance_pi(Ph);
      if (acc > 0){
	sidis.CalculateRfactor();
	if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	xQ2->Fill(x, Q2, weight * acc);
	xW->Fill(x, W, weight * acc);
	xz->Fill(x, z, weight * acc);
	xPt->Fill(x, Pt, weight * acc);
	xWp->Fill(x, Wp, weight * acc);
	zPt->Fill(z, Pt, weight * acc);
	zQ2->Fill(z, Q2, weight * acc);
	zW->Fill(z, W, weight * acc);
	zWp->Fill(z, Wp, weight * acc);
	PtQ2->Fill(Pt, Q2, weight * acc);
	PtW->Fill(Pt, W, weight * acc);
	PtWp->Fill(Pt, Wp, weight * acc);
	WQ2->Fill(W, Q2, weight * acc);
	WpQ2->Fill(Wp, Q2, weight * acc);
      }
    }
  }
  xQ2->Scale(lumi/Nsim);
  xW->Scale(lumi/Nsim);
  xz->Scale(lumi/Nsim);
  xPt->Scale(lumi/Nsim);
  xWp->Scale(lumi/Nsim);
  zW->Scale(lumi/Nsim);
  zQ2->Scale(lumi/Nsim);
  zWp->Scale(lumi/Nsim);
  zPt->Scale(lumi/Nsim);
  PtQ2->Scale(lumi/Nsim);
  WQ2->Scale(lumi/Nsim);
  PtW->Scale(lumi/Nsim);
  PtWp->Scale(lumi/Nsim);
  WpQ2->Scale(lumi/Nsim);
  fs->Write();
  return 0;
}

int MakeRateDistributionPlotZ(const double Ebeam){//Make z-? plot
  double lumi = 1.0e+10 * pow(0.197327, 2);
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2,1);
  sidis.SetHadron("pi+");
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double Xmin[6] = {0.0, 1.0, 0.01, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 10.0, 0.99, 2.0, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  gStyle->SetOptStat(0);
  //(Ph, z)
  TH2D * Phz = new TH2D("Phz", "", 500, 0.0, 10.0, 500, 0.0, 1.0);
  Phz->GetXaxis()->SetTitle("P_{h} / GeV");
  Phz->GetXaxis()->SetTitleSize(0.05);
  Phz->GetXaxis()->CenterTitle(true);
  Phz->GetXaxis()->SetLabelSize(0.05);
  Phz->GetYaxis()->SetTitle("z");
  Phz->GetYaxis()->SetTitleSize(0.05);
  Phz->GetYaxis()->CenterTitle(true);
  Phz->GetYaxis()->SetLabelSize(0.05);
  double x, Q2, z, Pt, W, Wp;
  double weight, acc;
  TLorentzVector lp, Ph;
  Long64_t Nsim = 100000000;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i % 10000000 == 0) std::cout << i * 100 / Nsim << " %" << std::endl;
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      z = sidis.GetVariable("z");
      //if (z < 0.3 || z > 0.7) continue;
      Q2 = sidis.GetVariable("Q2");
      if (Q2 < 1.0) continue;
      W = sidis.GetVariable("W");
      if (W < 2.3) continue;
      Wp = sidis.GetVariable("Wp");
      if (Wp < 1.6) continue;
      x = sidis.GetVariable("x");
      Pt = sidis.GetVariable("Pt");
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp, "all") * GetAcceptance_pi(Ph);
      if (acc > 0){
	//sidis.CalculateRfactor();
	//if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	Phz->Fill(Ph.P(), z, weight * acc);
      }
    }
  }
  Phz->Scale(lumi/Nsim);
  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  c0->SetLogz();
  Phz->DrawClone("colz");
  c0->Print("Phz.pdf");
  return 0;
}

int GenerateBinInfoFile(const char * filename, const double Ebeam, const char * hadron){//Bin the data and create the bin info file
  FILE * fp = fopen(filename, "w");
  fprintf(fp, "Q2l\t Q2u\t zl\t zu\t Ptl\t Ptu\t xl\t xu\n");
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2, 1);
  sidis.SetHadron(hadron);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+10 * pow(0.197327, 2);
  double eff = 0.85;
  double time = 48.0 * 24.0 * 3600.0;
  if (Ebeam < 10.0) time = 21.0 * 24.0 * 3600.0;
  double Nsim = 1.0e6;
  double Xmin[6] = {0.0, 0.0, 0.0, 0.0, -M_PI, -M_PI}; 
  double Xmax[6] = {0.7, 0.0, 0.0, 0.0, M_PI, M_PI};;//x, Q2, z, Pt, phih, phiS
  double weight = 0;
  double acc = 0;
  int Nx = 0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  double Q2list[7] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0};
  double statlist[6] = {1.9e7, 1.1e7, 5.0e6, 3.0e6, 2.0e6, 2.0e6};
  double zlist[9] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7};
  double Ptlist[7] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.6};
  int xi = 1;
  for (int Qi = 0; Qi < 6; Qi++){//Q2 loop
    Xmin[1] = Q2list[Qi];
    Xmax[1] = Q2list[Qi+1];
    for (int zi = 0; zi < 8; zi++){//z loop
      Xmin[2] = zlist[zi];
      Xmax[2] = zlist[zi+1];
      Xmin[3] = Ptlist[0];
      for (int kj = 1; kj < 7;){//Pt loop
	Xmax[3] = Ptlist[kj];
	sidis.SetRange(Xmin, Xmax);
	TH1D * hx = new TH1D("hx", "hx", 7000, 0.0, 0.7);
	printf("Q2:%.1f-%.1f  z:%.2f-%.2f  Pt:%.1f-%.1f\n",
	       Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3]);
	for (Long64_t i = 0; i < Nsim; i++){//generate events
	  weight = sidis.GenerateEvent(0, 1);
	  if (weight > 0){
            if (sidis.GetVariable("W") < 2.3) continue;
            if (sidis.GetVariable("Wp") < 1.6) continue;
	    sidis.CalculateRfactor();
	    if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	    lp = sidis.GetLorentzVector("lp");
	    Ph = sidis.GetLorentzVector("Ph");
	    acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph);
	    if (acc > 0)
	      hx->Fill(sidis.GetVariable("x"), weight * acc);
	  }
	}
	hx->Scale(lumi * time * eff / Nsim);
	if ((hx->Integral(1, -1) < statlist[Qi] && kj < 6) || (hx->Integral(1, -1) < 0.25 * statlist[Qi] && kj == 6)){
	  hx->Delete();
	  kj++;
	  continue;
	}
	Nx = 0;
	xi = 1;
	for (int xj = 1; xj <= 7000; xj++){
	  if (hx->Integral(xi, xj) > statlist[Qi] || xj == 7000){
	    fprintf(fp, "%.1f\t %.1f\t %.2f\t %.2f\t %.1f\t %.1f\t %.4f\t %.4f\n",
		    Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3],
		    hx->GetBinLowEdge(xi), hx->GetBinLowEdge(xj+1));
	    Nx++;
	    xi = xj + 1;
	  }
	}
	std::cout << Nx << std::endl;
	hx->Delete();
	Xmin[3] = Ptlist[kj];
	kj++;
      }
    }
  }
  fclose(fp);
  return 0;
}

int AnalyzeEstatUT3(const char * readfile, const char * savefile, const double Ebeam, const char * had){//bin analysis including stat. errors
  double Hadron = 0;
  if (strcmp(had, "pi+") == 0) Hadron = 0;
  else if (strcmp(had, "pi-") == 0) Hadron = 1;
  double Nucleon = 0;
  TFile * fs = new TFile(savefile, "RECREATE");
  TTree * Ts = new TTree("data", "data");
  Ts->SetDirectory(fs);
  double Eb = Ebeam;
  double x, y, z, Q2, Pt, phih, phiS;
  double dx, dy, dz, dQ2, dPt, dphih, dphiS, dv;
  double Nacc, fn;
  double Estatraw[3], Estat[3];
  Ts->Branch("Nucleon", &Nucleon, "Nucleon/D");
  Ts->Branch("Hadron", &Hadron, "Hadron/D");
  Ts->Branch("Ebeam", &Eb, "Ebeam/D");
  Ts->Branch("x", &x, "x/D");
  Ts->Branch("y", &y, "y/D");
  Ts->Branch("z", &z, "z/D");
  Ts->Branch("Q2", &Q2, "Q2/D");
  Ts->Branch("Pt", &Pt, "Pt/D");
  Ts->Branch("dx", &dx, "dx/D");
  Ts->Branch("dy", &dy, "dy/D");
  Ts->Branch("dz", &dz, "dz/D");
  Ts->Branch("dQ2", &dQ2, "dQ2/D");
  Ts->Branch("dPt", &dPt, "dPt/D");
  Ts->Branch("dphih", &dphih, "dphih/D");
  Ts->Branch("dphiS", &dphiS, "dphiS/D");
  Ts->Branch("dv", &dv, "dv/D");
  Ts->Branch("Nacc", &Nacc, "Nacc/D");
  Ts->Branch("fn", &fn, "fn/D");
  Ts->Branch("E0statraw", &Estatraw[0], "E0statraw/D");
  Ts->Branch("E1statraw", &Estatraw[1], "E1statraw/D");
  Ts->Branch("E2statraw", &Estatraw[2], "E2statraw/D");
  Ts->Branch("E0stat", &Estat[0], "E0stat/D");
  Ts->Branch("E1stat", &Estat[1], "E1stat/D");
  Ts->Branch("E2stat", &Estat[2], "E2stat/D");
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2, 1);
  sidis.SetHadron(had);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+10 * pow(0.197327, 2);
  double eff = 0.85;
  double time = 48.0 * 24.0 * 3600.0;
  if (Ebeam < 10.0) time = 21.0 * 24.0 * 3600.0;
  Long64_t Nsim = 0;
  Long64_t Nrec = 0;
  double Xmin[6] = {0.0, 0.0, 0.0, 0.0, -M_PI, -M_PI}; 
  double Xmax[6] = {0.7, 0.0, 0.0, 0.0, M_PI, M_PI};;//x, Q2, z, Pt, phih, phiS
  double weight = 0;
  double weight_n = 0;
  double acc = 0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);
  Lsidis sidis_n;
  sidis_n.SetNucleus(0, 1);
  sidis_n.SetHadron(had);
  sidis_n.SetInitialState(l, P);
  sidis_n.SetPDFset("CT14lo");
  sidis_n.SetFFset("DSSFFlo");
  ifstream infile(readfile);
  char tmp[300];
  infile.getline(tmp, 256);
  int Nt = 0;
  while (infile >> Xmin[1] >> Xmax[1] >> Xmin[2] >> Xmax[2] >> Xmin[3] >> Xmax[3] >> Xmin[0] >> Xmax[0]){
    printf("%.4d  Q2[%.1f,%.1f]  z[%.2f,%.2f]  Pt[%.1f,%.1f]  x[%.4f,%.4f]\n",
	   Nt++, Xmin[1], Xmax[1], Xmin[2], Xmax[2], Xmin[3], Xmax[3], Xmin[0], Xmax[0]);
    sidis.SetRange(Xmin, Xmax);
    sidis_n.SetRange(Xmin, Xmax);
    TH1D * hvar = new TH1D("hvar", "hvar", 7, -0.5, 6.5);
    TH2D * hs = new TH2D("hs", "hs", 36, -M_PI, M_PI, 18, 0, M_PI);
    Nsim = 0;
    Nrec = 0;
    for (Long64_t i = 0; i < 1.0e7; i++){
      Nsim++;
      weight = sidis.GenerateEvent(0, 1);
      if (weight > 0){
        if (sidis.GetVariable("W") < 2.3) continue;
        if (sidis.GetVariable("Wp") < 1.6) continue;
	sidis.CalculateRfactor();
	if (sidis.GetVariable("Rfactor") > Rfactor0) continue;
	lp = sidis.GetLorentzVector("lp");
	Ph = sidis.GetLorentzVector("Ph");
	acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph);
	if (acc > 0){
	  sidis_n.SetFinalState(lp, Ph);
	  sidis_n.CalculateVariables();
	  weight_n = sidis_n.GetEventWeight(0, 1);
	  Nrec++;
	  hvar->Fill(0., weight_n * acc);
	  hvar->Fill(1., weight * acc);
	  hvar->Fill(2., weight * acc * sidis.GetVariable("x"));
	  hvar->Fill(3., weight * acc * sidis.GetVariable("y"));
	  hvar->Fill(4., weight * acc * sidis.GetVariable("z"));
	  hvar->Fill(5., weight * acc * sidis.GetVariable("Q2"));
	  hvar->Fill(6., weight * acc * sidis.GetVariable("Pt"));
	  hs->Fill(sidis.GetVariable("phih"), std::abs(sidis.GetVariable("phiS")), weight * acc);
	}
      }
      if (Nrec > 500000) break;
    }
    hvar->Scale(lumi * time * eff / Nsim);
    hs->Scale(lumi * time * eff / Nsim);
    Nacc = hvar->GetBinContent(2);
    fn = hvar->GetBinContent(1) / Nacc;
    x = hvar->GetBinContent(3) / Nacc;
    y = hvar->GetBinContent(4) / Nacc;
    z = hvar->GetBinContent(5) / Nacc;
    Q2 = hvar->GetBinContent(6) / Nacc;
    Pt = hvar->GetBinContent(7) / Nacc;
    TMatrixD MUT3(3,3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	MUT3(i,j) = 0.0;
    for (int i = 1; i <= 36; i++){
      for (int j = 1; j <= 18; j++){
	phih = hs->GetXaxis()->GetBinCenter(i);
	phiS = hs->GetYaxis()->GetBinCenter(j);
	MUT3(0,0) += sin(phih - phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(0,1) += sin(phih - phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(0,2) += sin(phih - phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,0) += sin(phih + phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,1) += sin(phih + phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(1,2) += sin(phih + phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,0) += sin(3.0 * phih - phiS) * sin(phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,1) += sin(3.0 * phih - phiS) * sin(phih + phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
	MUT3(2,2) += sin(3.0 * phih - phiS) * sin(3.0 * phih - phiS) * hs->GetBinContent(i,j) / Nacc * 2.0 * M_PI * M_PI;
      }
    }
    MUT3.Invert();
    for (int i = 0; i < 3; i++){
      Estatraw[i] = sqrt(2.0 * M_PI * M_PI / Nacc * (pow(MUT3(i,0),2) + pow(MUT3(i,1), 2) + pow(MUT3(i,2), 2)) * M_PI * M_PI);
      Estat[i] = Estatraw[i] / fn / 0.6 / 0.86;
      if (isnan(Estat[i]))
	std::cout << "NaN warning in Estat!" << std::endl;
    }
    Ts->Fill();
    hvar->Delete();
    hs->Delete();
  }
  fs->Write();
  infile.close();
  return 0;
}
  
double CheckCurrentCut(const double Ebeam, const char * hadron, const double kT2 = 0.16, const double MiT2 = 0.4, const double MfT2 = 0.4, const char * plotname = 0){
  Lsidis sidis;
  TLorentzVector l(0, 0, Ebeam, Ebeam);
  TLorentzVector P(0, 0, 0, 0.938272);
  sidis.SetNucleus(2, 1);
  sidis.SetHadron(hadron);
  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CT14lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+10 * pow(0.197327, 2);
  double Nsim = 1.0e7;
  TH2D * h0 = new TH2D("h0", "", 1, 0.2, 0.8, 1, 0.0, 1.6);
  h0->GetXaxis()->SetTitle("z");
  h0->GetXaxis()->CenterTitle(true);
  h0->GetXaxis()->SetTitleSize(0.05);
  h0->GetXaxis()->SetTitleOffset(1.15);
  h0->GetXaxis()->SetLabelSize(0.055);
  h0->GetYaxis()->SetTitle("P_{hT} / GeV");
  h0->GetYaxis()->CenterTitle(true);
  h0->GetYaxis()->SetTitleSize(0.05);
  h0->GetYaxis()->SetTitleOffset(1.15);
  h0->GetYaxis()->SetLabelSize(0.055);
  TH2D * hall = new TH2D("hall", "Before cut", 60, 0.2, 0.8, 160, 0.0, 1.6);
  TH2D * hcut = new TH2D("hcut", "After cut", 60, 0.2, 0.8, 160, 0.0, 1.6);
  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 8.0, 0.7, 1.6, M_PI, M_PI};
  sidis.SetRange(Xmin, Xmax);
  double weight = 0;
  double acc = 0;
  TLorentzVector lp, Ph;
  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 1);
    if (weight > 0){
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph);
      if (acc > 0){
	hall->Fill(sidis.GetVariable("z"), sidis.GetVariable("Pt"), weight * acc);
	sidis.CalculateRfactor(kT2, MiT2, MfT2);
	if (sidis.GetVariable("Rfactor") < 0.4){
	  hcut->Fill(sidis.GetVariable("z"), sidis.GetVariable("Pt"), weight * acc);
	}
      }
    }
  }
  hall->Scale(lumi/Nsim);
  hcut->Scale(lumi/Nsim);
  double rate = hcut->Integral(1, -1);
  std::cout << "All: " << hall->Integral(1, -1) << "   Cut: " << hcut->Integral(1, -1) << std::endl;
  if (plotname != 0){
    gStyle->SetOptStat(0);
    //hall->GetZaxis()->SetRangeUser(0.01, hall->GetMaximum()/0.95);
    hcut->GetZaxis()->SetRangeUser(0.01, hall->GetMaximum());
    TCanvas * c0 = new TCanvas("c0", "", 1600, 600);
    c0->SetBorderMode(0);
    c0->SetBorderSize(2);
    c0->SetFrameBorderMode(0);
    c0->Divide(2, 1);
    c0->cd(1);
    c0->cd(1)->SetLeftMargin(0.15);
    c0->cd(1)->SetBottomMargin(0.15);
    h0->Draw();
    hall->Draw("samecolz");
    c0->cd(2);
    c0->cd(2)->SetLeftMargin(0.15);
    c0->cd(2)->SetBottomMargin(0.15);
    h0->Draw();
    hcut->Draw("samecolz");
    c0->Print(plotname);
    c0->Close();
  }
  h0->Delete();
  hall->Delete();
  hcut->Delete();
  return rate;
}

int CreateFileSivers(const char * rootfile1, const char * rootfile2, const char * csvfile){//Create file for Sivers analysis use
  TChain * Ts = new TChain("data", "data");
  Ts->Add(rootfile1);
  Ts->Add(rootfile2);
  double Nucleon, Hadron, Ebeam, x, y, z, Q2, Pt, stat, systrel, systabs, fn;
  Ts->SetBranchAddress("Nucleon", &Nucleon);
  Ts->SetBranchAddress("Hadron", &Hadron);
  Ts->SetBranchAddress("Ebeam", &Ebeam);
  Ts->SetBranchAddress("x", &x);
  Ts->SetBranchAddress("y", &y);
  Ts->SetBranchAddress("z", &z);
  Ts->SetBranchAddress("Q2", &Q2);
  Ts->SetBranchAddress("Pt", &Pt);
  Ts->SetBranchAddress("E1stat", &stat);
  Ts->SetBranchAddress("fn", &fn);
  FILE * file = fopen(csvfile, "w");
  fprintf(file, "i,Ebeam,x,y,z,Q2,pT,obs,value,stat,systrel,systabs,target,hadron,Experiment\n");
  for (int i = 0; i < Ts->GetEntries(); i++){
    std::cout << i << std::endl;
    Ts->GetEntry(i);
    systrel = 0.0;
    systabs = 0.0;
    systrel += pow(0.03, 2);//target polarization
    systrel += pow(0.05, 2);//nuclear effect
    systrel += pow(0.025, 2);//radiative correction
    systrel += pow(0.03, 2);//diffractive meson
    systrel += pow(0.002, 2);//random coincidence
    if (Ebeam > 10.0)//raw asymmetry
      systabs += 1.7e-4 / 0.6 / fn / 0.86;
    else
      systabs += 2.57e-4 / 0.6 / fn / 0.86;
    systrel = sqrt(systrel);
    if (Hadron == 0)
      fprintf(file, "%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      i, Ebeam, x, y, z, Q2, Pt, "AUT", 0.0, stat, systrel, systabs, "neutron", "pi+", "solid");
    else if (Hadron == 1)
      fprintf(file, "%d,%.1f,%.6f,%.6f,%.6f,%.6f,%.6f,%s,%.1f,%.6f,%.6f,%.6f,%s,%s,%s\n",
	      i, Ebeam, x, y, z, Q2, Pt, "AUT", 0.0, stat, systrel, systabs, "neutron", "pi-", "solid");
  }
  fclose(file);
  return 0;
}







  





  
	
	
