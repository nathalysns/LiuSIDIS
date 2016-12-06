#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../Header/Lsidis1.h"

using namespace std;

TFile * file_negative = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root","r");
TFile * file_positive = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_positive_output.root","r");
TH2F * acc_FA_negative = (TH2F *) file_negative->Get("acceptance_forwardangle");
TH2F * acc_LA_negative = (TH2F *) file_negative->Get("acceptance_largeangle");
TH2F * acc_FA_positive = (TH2F *) file_positive->Get("acceptance_forwardangle");
TH2F * acc_LA_positive = (TH2F *) file_positive->Get("acceptance_largeangle");

double GetAcceptance(const TLorentzVector p, const char * part){//
  double theta = p.Theta() / M_PI * 180.0;
  if (theta < 8.0) return 0;
  double mom = p.P();
  double acc = 0;
  if (strcmp(part, "e-") == 0){
    acc += acc_FA_negative->GetBinContent(acc_FA_negative->GetXaxis()->FindBin(theta), acc_FA_negative->GetYaxis()->FindBin(mom));
    if (mom > 3.5){
      acc += acc_LA_negative->GetBinContent(acc_LA_negative->GetXaxis()->FindBin(theta), acc_LA_negative->GetYaxis()->FindBin(mom));
    }
  }
  else if (strcmp(part, "pi+") == 0){
    acc += acc_FA_positive->GetBinContent(acc_FA_positive->GetXaxis()->FindBin(theta), acc_FA_positive->GetYaxis()->FindBin(mom));
  }
  else if (strcmp(part, "pi-") == 0){
    acc += acc_FA_negative->GetBinContent(acc_FA_negative->GetXaxis()->FindBin(theta), acc_FA_negative->GetYaxis()->FindBin(mom));
  }
  return acc;
}

int main(int argc, char * argv[]){

  if (argc < 2){
    cout << "missing hadron input" << endl;
    return -1;
  };
  
  gRandom->SetSeed(0);

  Lsidis mysidis;

  TLorentzVector l(0, 0, 11.0, 11.0);
  TLorentzVector P(0, 0, 0, 0.938272);

  double lumi = 1.0e+10 * pow(0.197327, 2);
  mysidis.SetNucleus(2, 1);
  mysidis.SetHadron(argv[1]);
  
  mysidis.SetInitialState(l, P);
  mysidis.SetPDFset("CT14lo");
  mysidis.SetFFset("DSSFFlo");

  double Xmin[6] = {0.03, 1.0, 0.3, 0.0, -M_PI, -M_PI};//x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.75, 10.0, 0.7, 1.8, M_PI, M_PI};
  mysidis.SetRange(Xmin, Xmax);

  double weight = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);

  double Nsim = 1.0e8;

  TH1D * tx = new TH1D("tx", "", 70, 0.0, 0.7);
  TH1D * sx = new TH1D("sx", "", 70, 0.0, 0.7);
  TH1D * rx = new TH1D("rx", "", 70, 0.0, 0.7);

  TH1D * tz = new TH1D("tz", "", 50, 0.25, 0.75);
  TH1D * sz = new TH1D("sz", "", 50, 0.25, 0.75);
  TH1D * rz = new TH1D("rz", "", 50, 0.25, 0.75);

  TH1D * tPt = new TH1D("tPt", "", 60, 0.0, 1.8);
  TH1D * sPt = new TH1D("sPt", "", 60, 0.0, 1.8);
  TH1D * rPt = new TH1D("rPt", "", 60, 0.0, 1.8);

  TH2D * txz = new TH2D("txz", "", 35, 0.0, 0.7, 30, 0.25, 0.75);
  TH2D * sxz = new TH2D("sxz", "", 35, 0.0, 0.7, 30, 0.25, 0.75);
  TH2D * rxz = new TH2D("rxz", "", 35, 0.0, 0.7, 30, 0.25, 0.75);

  TH2D * txPt = new TH2D("txPt", "", 35, 0.0, 0.7, 36, 0.0, 1.8);
  TH2D * sxPt = new TH2D("sxPt", "", 35, 0.0, 0.7, 36, 0.0, 1.8);
  TH2D * rxPt = new TH2D("rxPt", "", 35, 0.0, 0.7, 36, 0.0, 1.8);

  TH2D * tzPt = new TH2D("tzPt", "", 30, 0.25, 0.75, 36, 0.0, 1.8);
  TH2D * szPt = new TH2D("szPt", "", 30, 0.25, 0.75, 36, 0.0, 1.8);
  TH2D * rzPt = new TH2D("rzPt", "", 30, 0.25, 0.75, 36, 0.0, 1.8);

  double x, Q2, z, Pt, W, Wp;
  double acc, th;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%10000000 == 0) cout << i << endl;
    weight = mysidis.GenerateEvent(0, 1);
    if (weight > 0){
      lp = mysidis.GetLorentzVector("lp");
      Ph = mysidis.GetLorentzVector("Ph");

      x = mysidis.GetVariable("x");
      Q2 = mysidis.GetVariable("Q2");
      z = mysidis.GetVariable("z");
      Pt = mysidis.GetVariable("Pt");
      W = mysidis.GetVariable("W");
      Wp = mysidis.GetVariable("Wp");

      if (Q2 < 1.0) continue;
      if (z < 0.3 || z > 0.7) continue;
      
      if (W < 2.3) continue;
      if (Wp < 1.6) continue;

      acc = GetAcceptance(lp, "e-") * GetAcceptance(Ph, argv[1]);
      if (acc == 0) continue;

      th = lp.Theta() * 180.0 / M_PI;

      tx->Fill(x, weight * acc);
      tz->Fill(z, weight * acc);
      tPt->Fill(Pt, weight * acc);
      txz->Fill(x, z, weight * acc);
      txPt->Fill(x, Pt, weight * acc);
      tzPt->Fill(z, Pt, weight * acc);

      if ((th >= 8.0 && th < 9.0 && lp.P() < 5.0) || (th >= 9.0 && th < 10.0 && lp.P() < 4.0) || (th >= 10.0 && th < 12.0 && lp.P() < 3.0) || (th >= 12.0 && th < 14.0 && lp.P() < 2.0) || (th >= 14.0 && th < 15.0 && lp.P() < 1.0) ){
	sx->Fill(x, weight * acc / 10.0);
	sz->Fill(z, weight * acc / 10.0);
	sPt->Fill(Pt, weight * acc / 10.0);
	sxz->Fill(x, z, weight * acc / 10.0);
	sxPt->Fill(x, Pt, weight * acc / 10.0);
	szPt->Fill(z, Pt, weight * acc / 10.0);
	rx->Fill(x, weight * acc);
	rz->Fill(z, weight * acc);
	rPt->Fill(Pt, weight * acc);
	rxz->Fill(x, z, weight * acc);
	rxPt->Fill(x, Pt, weight * acc);
	rzPt->Fill(z, Pt, weight * acc);
      }
      else{
	sx->Fill(x, weight * acc);
	sz->Fill(z, weight * acc);
	sPt->Fill(Pt, weight * acc);
	sxz->Fill(x, z, weight * acc);
	sxPt->Fill(x, Pt, weight * acc);
	szPt->Fill(z, Pt, weight * acc);
      }
    }
  }

  tx->Scale(lumi/Nsim);
  tz->Scale(lumi/Nsim);
  tPt->Scale(lumi/Nsim);
  txz->Scale(lumi/Nsim);
  txPt->Scale(lumi/Nsim);
  tzPt->Scale(lumi/Nsim);
  sx->Scale(lumi/Nsim);
  sz->Scale(lumi/Nsim);
  sPt->Scale(lumi/Nsim);
  sxz->Scale(lumi/Nsim);
  sxPt->Scale(lumi/Nsim);
  szPt->Scale(lumi/Nsim);
  rx->Scale(lumi/Nsim);
  rz->Scale(lumi/Nsim);
  rPt->Scale(lumi/Nsim);
  rxz->Scale(lumi/Nsim);
  rxPt->Scale(lumi/Nsim);
  rzPt->Scale(lumi/Nsim);

  ////
  sx->SetLineColor(2);
  sx->SetLineWidth(1.5);
  sx->GetXaxis()->SetTitle("x");
  sx->GetXaxis()->SetTitleSize(0.05);
  sx->GetXaxis()->SetLabelSize(0.05);
  sx->GetYaxis()->SetLabelSize(0.05);
  tx->SetLineColor(4);
  tx->SetLineWidth(1.5);
  tx->GetXaxis()->SetTitle("x");
  tx->GetXaxis()->SetTitleSize(0.05);
  tx->GetXaxis()->SetLabelSize(0.05);
  tx->GetYaxis()->SetLabelSize(0.05);

  sz->SetLineColor(2);
  sz->SetLineWidth(1.5);
  sz->GetXaxis()->SetTitle("z");
  sz->GetXaxis()->SetTitleSize(0.05);
  sz->GetXaxis()->SetLabelSize(0.05);
  sz->GetYaxis()->SetLabelSize(0.05);
  tz->SetLineColor(4);
  tz->SetLineWidth(1.5);
  tz->GetXaxis()->SetTitle("z");
  tz->GetXaxis()->SetTitleSize(0.05);
  tz->GetXaxis()->SetLabelSize(0.05);
  tz->GetYaxis()->SetLabelSize(0.05);

  sPt->SetLineColor(2);
  sPt->SetLineWidth(1.5);
  sPt->GetXaxis()->SetTitle("P_{hT} (GeV)");
  sPt->GetXaxis()->SetTitleSize(0.05);
  sPt->GetXaxis()->SetLabelSize(0.05);
  sPt->GetYaxis()->SetLabelSize(0.05);
  tPt->SetLineColor(4);
  tPt->SetLineWidth(1.5);
  tPt->GetXaxis()->SetTitle("P_{hT} (GeV)");
  tPt->GetXaxis()->SetTitleSize(0.05);
  tPt->GetXaxis()->SetLabelSize(0.05);
  tPt->GetYaxis()->SetLabelSize(0.05);

  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1", "", 1200, 1800);
  c1->Divide(2, 3);
  c1->cd(1);
  tx->DrawClone();
  sx->DrawClone("same");
  c1->cd(2);
  sx->Divide(tx);
  sx->GetYaxis()->SetRangeUser(0, 1.2);
  sx->DrawClone();
  c1->cd(3);
  tz->DrawClone();
  sz->DrawClone("same");
  c1->cd(4);
  sz->Divide(tz);
  sz->GetYaxis()->SetRangeUser(0, 1.2);
  sz->DrawClone();
  c1->cd(5);
  tPt->DrawClone();
  sPt->DrawClone("same");
  c1->cd(6);
  sPt->Divide(tPt);
  sPt->GetYaxis()->SetRangeUser(0, 1.2);
  sPt->DrawClone();

  c1->Print("Prescale10(P5433221)1D.pdf");

  ////
  sxz->GetXaxis()->SetTitle("x");
  sxz->GetYaxis()->SetTitle("z");
  sxz->GetXaxis()->SetTitleSize(0.05);
  sxz->GetYaxis()->SetTitleSize(0.05);
  sxz->GetXaxis()->SetLabelSize(0.05);
  sxz->GetYaxis()->SetLabelSize(0.05);

  sxPt->GetXaxis()->SetTitle("x");
  sxPt->GetYaxis()->SetTitle("P_{hT} (GeV)");
  sxPt->GetXaxis()->SetTitleSize(0.05);
  sxPt->GetYaxis()->SetTitleSize(0.05);
  sxPt->GetXaxis()->SetLabelSize(0.05);
  sxPt->GetYaxis()->SetLabelSize(0.05);

  szPt->GetXaxis()->SetTitle("z");
  szPt->GetYaxis()->SetTitle("P_{hT} (GeV)");
  szPt->GetXaxis()->SetTitleSize(0.05);
  szPt->GetYaxis()->SetTitleSize(0.05);
  szPt->GetXaxis()->SetLabelSize(0.05);
  szPt->GetYaxis()->SetLabelSize(0.05);

  TCanvas * c2 = new TCanvas("c2", "", 600, 1800);
  c2->Divide(1, 3);
  c2->cd(1);
  sxz->Divide(txz);
  sxz->DrawClone("colz");
  c2->cd(2);
  sxPt->Divide(txPt);
  sxPt->DrawClone("colz");
  c2->cd(3);
  szPt->Divide(tzPt);
  szPt->DrawClone("colz");

  c2->Print("Prescale10(P5433221)2D.pdf");

  ////
  rxz->GetXaxis()->SetTitle("x");
  rxz->GetYaxis()->SetTitle("z");
  rxz->GetXaxis()->SetTitleSize(0.05);
  rxz->GetYaxis()->SetTitleSize(0.05);
  rxz->GetXaxis()->SetLabelSize(0.05);
  rxz->GetYaxis()->SetLabelSize(0.05);
  
  rxPt->GetXaxis()->SetTitle("x");
  rxPt->GetYaxis()->SetTitle("P_{hT} (GeV)");
  rxPt->GetXaxis()->SetTitleSize(0.05);
  rxPt->GetYaxis()->SetTitleSize(0.05);
  rxPt->GetXaxis()->SetLabelSize(0.05);
  rxPt->GetYaxis()->SetLabelSize(0.05);

  rzPt->GetXaxis()->SetTitle("z");
  rzPt->GetYaxis()->SetTitle("P_{hT} (GeV)");
  rzPt->GetXaxis()->SetTitleSize(0.05);
  rzPt->GetYaxis()->SetTitleSize(0.05);
  rzPt->GetXaxis()->SetLabelSize(0.05);
  rzPt->GetYaxis()->SetLabelSize(0.05);

  TCanvas * c3 = new TCanvas("c3", "", 600, 1800);
  c3->Divide(1, 3);
  c3->cd(1);
  rxz->DrawClone("colz");
  c3->cd(2);
  rxPt->DrawClone("colz");
  c3->cd(3);
  rzPt->DrawClone("colz");

  c3->Print("Distribution(P5433221)2D.pdf");

  cout << tx->Integral(1, -1) << endl;


  return 0;
}
