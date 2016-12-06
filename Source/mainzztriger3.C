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

  TH2D * t8 = new TH2D("t8", "", 70, 0.0, 0.7, 80, 0.0, 8.0);
  TH2D * t9 = new TH2D("t9", "", 70, 0.0, 0.7, 80, 0.0, 8.0);
  TH2D * t10 = new TH2D("t10", "", 70, 0.0, 0.7, 80, 0.0, 8.0);
  TH2D * t11 = new TH2D("t11", "", 70, 0.0, 0.7, 80, 0.0, 8.0);
  TH2D * t12 = new TH2D("t12", "", 70, 0.0, 0.7, 80, 0.0, 8.0);
  TH2D * t13 = new TH2D("t13", "", 70, 0.0, 0.7, 80, 0.0, 8.0);

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
      if (th >= 8.0 && th < 9.0)
	t8->Fill(x, lp.P(), weight * acc);
      else if (th >= 9.0 && th < 10.0)
	t9->Fill(x, lp.P(), weight * acc);
      else if (th >= 10.0 && th < 11.0)
	t10->Fill(x, lp.P(), weight * acc);
      else if (th >= 11.0 && th < 12.0)
	t11->Fill(x, lp.P(), weight * acc);
      else if (th >= 12.0 && th < 13.0)
	t12->Fill(x, lp.P(), weight * acc);
      else if (th >= 13.0 && th < 14.0)
	t13->Fill(x, lp.P(), weight * acc);
    }
  }

  t8->Scale(lumi/Nsim);
  t9->Scale(lumi/Nsim);
  t10->Scale(lumi/Nsim);
  t11->Scale(lumi/Nsim);
  t12->Scale(lumi/Nsim);
  t13->Scale(lumi/Nsim);
 
   ////
  t8->SetTitle("Theta: 8-9 Deg");
  t8->GetXaxis()->SetTitle("x");
  t8->GetYaxis()->SetTitle("P_{e} (GeV)");
  t8->GetXaxis()->SetTitleSize(0.05);
  t8->GetYaxis()->SetTitleSize(0.05);
  t8->GetXaxis()->SetLabelSize(0.05);
  t8->GetYaxis()->SetLabelSize(0.05);

  t9->SetTitle("Theta: 9-10 Deg");
  t9->GetXaxis()->SetTitle("x");
  t9->GetYaxis()->SetTitle("P_{e} (GeV)");  
  t9->GetXaxis()->SetTitleSize(0.05);                                      
  t9->GetYaxis()->SetTitleSize(0.05);
  t9->GetXaxis()->SetLabelSize(0.05);
  t9->GetYaxis()->SetLabelSize(0.05);

  t10->SetTitle("Theta: 10-11 Deg");
  t10->GetXaxis()->SetTitle("x");
  t10->GetYaxis()->SetTitle("P_{e} (GeV)");
  t10->GetXaxis()->SetTitleSize(0.05);
  t10->GetYaxis()->SetTitleSize(0.05);
  t10->GetXaxis()->SetLabelSize(0.05);
  t10->GetYaxis()->SetLabelSize(0.05);

  t11->SetTitle("Theta: 11-12 Deg");
  t11->GetXaxis()->SetTitle("x");
  t11->GetYaxis()->SetTitle("P_{e} (GeV)");
  t11->GetXaxis()->SetTitleSize(0.05);
  t11->GetYaxis()->SetTitleSize(0.05);
  t11->GetXaxis()->SetLabelSize(0.05);
  t11->GetYaxis()->SetLabelSize(0.05);

  t12->SetTitle("Theta: 12-13 Deg");
  t12->GetXaxis()->SetTitle("x");
  t12->GetYaxis()->SetTitle("P_{e} (GeV)");
  t12->GetXaxis()->SetTitleSize(0.05);
  t12->GetYaxis()->SetTitleSize(0.05);
  t12->GetXaxis()->SetLabelSize(0.05);
  t12->GetYaxis()->SetLabelSize(0.05);

  t13->SetTitle("Theta: 13-14 Deg");
  t13->GetXaxis()->SetTitle("x");
  t13->GetYaxis()->SetTitle("P_{e} (GeV)");
  t13->GetXaxis()->SetTitleSize(0.05);
  t13->GetYaxis()->SetTitleSize(0.05);
  t13->GetXaxis()->SetLabelSize(0.05);
  t13->GetYaxis()->SetLabelSize(0.05);

  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1", "", 1200, 1800);
  c1->Divide(2,3);
  c1->cd(1);
  t8->DrawClone("colz");
  c1->cd(2);
  t9->DrawClone("colz");
  c1->cd(3);
  t10->DrawClone("colz");
  c1->cd(4);
  t11->DrawClone("colz");
  c1->cd(5);
  t12->DrawClone("colz");
  c1->cd(6);
  t13->DrawClone("colz"); 

  c1->Print("Distribution(x,Pe)2D.pdf");

  return 0;
}
