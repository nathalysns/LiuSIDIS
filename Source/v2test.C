#include "Lsidis.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

TFile * file_negative = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root","r");
TFile * file_positive = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_positive_output.root","r");
TH2F * acc_FA_negative = (TH2F *) file_negative->Get("acceptance_forwardangle");
TH2F * acc_LA_negative = (TH2F *) file_negative->Get("acceptance_largeangle");
TH2F * acc_FA_positive = (TH2F *) file_positive->Get("acceptance_forwardangle");
TH2F * acc_LA_positive = (TH2F *) file_positive->Get("acceptance_largeangle");

double GetAcceptance(const TLorentzVector p, const char * dfile){//
  double theta = p.Theta() / M_PI * 180.0;
  double mom = p.P();
  double acc = 0;
  if (theta > 8.0 && mom > 0.8){
    if (strcmp(dfile, "e-") == 0){
      acc += acc_FA_negative->GetBinContent(acc_FA_negative->GetXaxis()->FindBin(theta), acc_FA_negative->GetYaxis()->FindBin(mom));
      if (mom > 3.5){
	acc += acc_LA_negative->GetBinContent(acc_LA_negative->GetXaxis()->FindBin(theta), acc_LA_negative->GetYaxis()->FindBin(mom));
      }
    }
    else if (strcmp(dfile, "h+") == 0){
      acc += acc_FA_positive->GetBinContent(acc_FA_positive->GetXaxis()->FindBin(theta), acc_FA_positive->GetYaxis()->FindBin(mom));
    }    
    else if (strcmp(dfile, "h-") == 0){
      acc += acc_FA_negative->GetBinContent(acc_FA_negative->GetXaxis()->FindBin(theta), acc_FA_negative->GetYaxis()->FindBin(mom));
    }
  }
  return acc;
}

int main(int argc, char * argv[]){

  if (argc < 2) {
    cout << "missing inputs" << endl;
    return -1;
  }

  gRandom->SetSeed(1);
  gStyle->SetOptStat(0);

  Lsidis mysidis;
  TLorentzVector l(0, 0, 11.0, 11.0);
  TLorentzVector P(0, 0, 0, 0.938272);

  mysidis.SetNucleus(2, 1);//helium-3
  mysidis.SetHadron(argv[1]);

  int ic = mysidis.GetHadronCharge();

  mysidis.SetInitialState(l, P);
  mysidis.SetPDFset("CT14lo");
  mysidis.SetFFset("DSSFFlo");

  double Q2min = 1.0;
  double Q2max = 10.0;
  double zmin = 0.3;
  double zmax = 0.7;
  
  double Xmin[6] = {Q2min/24.0, Q2min, zmin, 0.0, -M_PI, -M_PI};//x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.75, Q2max, zmax, 1.6, M_PI, M_PI};
  mysidis.SetRange(Xmin, Xmax);

  cout << "Q2: " << Q2min << " -- " << Q2max << "     z: " << zmin << " -- " << zmax << endl; 

  double weight = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);

  double lumi = 1.0e+10 * pow(0.197327, 2);

  
  Long64_t Nsim1 = 1000000;
  Long64_t Nsim2 = 10000;
  double acc = 0.0;

  TH1D * hx = new TH1D("hx", "", 100, 0.0, 0.7);
  TH1D * hQ2 = new TH1D("hQ2", "", 100, 0.0, 10.0);
  TH1D * hz = new TH1D("hz", "", 100, 0.25, 0.75);
  TH1D * hPt = new TH1D("hPt", "", 100, 0.0, 1.6);
  TH1D * hphih = new TH1D("hphih", "", 100, -3.2, 3.2);

  hx->GetXaxis()->SetTitle("x");
  hx->GetXaxis()->SetTitleSize(0.05);
  hx->GetXaxis()->SetLabelSize(0.05);
  hx->GetYaxis()->SetLabelSize(0.05);
  hx->GetXaxis()->CenterTitle();
  hQ2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  hQ2->GetXaxis()->SetTitleSize(0.05);
  hQ2->GetXaxis()->SetLabelSize(0.05);
  hQ2->GetYaxis()->SetLabelSize(0.05);
  hQ2->GetXaxis()->CenterTitle();
  hz->GetXaxis()->SetTitle("z");
  hz->GetXaxis()->SetTitleSize(0.05);
  hz->GetXaxis()->SetLabelSize(0.05);
  hz->GetYaxis()->SetLabelSize(0.05);
  hz->GetXaxis()->CenterTitle();
  hPt->GetXaxis()->SetTitle("P_{hT} (GeV)");
  hPt->GetXaxis()->SetTitleSize(0.05);
  hPt->GetXaxis()->SetLabelSize(0.05);
  hPt->GetYaxis()->SetLabelSize(0.05);
  hPt->GetXaxis()->CenterTitle();
  hphih->GetXaxis()->SetTitle("#phi_{h} (rad)");
  hphih->GetXaxis()->SetTitleSize(0.05);
  hphih->GetXaxis()->SetLabelSize(0.05);
  hphih->GetYaxis()->SetLabelSize(0.05);
  hphih->GetXaxis()->CenterTitle();

  hx->SetLineColor(4);
  hQ2->SetLineColor(4);
  hz->SetLineColor(4);
  hPt->SetLineColor(4);
  hphih->SetLineColor(4);
 
  double x, Q2, z, W, Wp, Pt, phih;

  TCanvas * c00 = new TCanvas("c00", "", 1200, 600);
  c00->Divide(3,2);
  for (Long64_t i = 0; i < Nsim1; i++){
    if (i%10000 == 0) cout << i << endl;
    weight = mysidis.GenerateEvent(0, 1);
    if (weight > 0){
      acc = 0;
      lp = mysidis.GetLorentzVector("lp");
      Ph = mysidis.GetLorentzVector("Ph");
      x = mysidis.GetVariable("x");
      Q2 = mysidis.GetVariable("Q2");
      z = mysidis.GetVariable("z");
      Pt = mysidis.GetVariable("Pt");
      W = mysidis.GetVariable("W");
      Wp = mysidis.GetVariable("Wp");
      phih = mysidis.GetVariable("phih");

      if (Q2 < 1.0) continue;
      if (W < 2.3) continue;
      if (Wp < 1.6) continue;

      if (ic == 1)
	acc = GetAcceptance(lp, "e-") * GetAcceptance(Ph, "h+");
      else if (ic == -1)
	acc = GetAcceptance(lp, "e-") * GetAcceptance(Ph, "h-");

      if (acc > 0){
	hx->Fill(x, weight * acc);
	hQ2->Fill(Q2, weight * acc);
	hz->Fill(z, weight * acc);
	hPt->Fill(Pt, weight * acc);
	hphih->Fill(phih, weight * acc);
      }
    }
  }

  hx->Scale(lumi/Nsim1);
  hQ2->Scale(lumi/Nsim1);
  hz->Scale(lumi/Nsim1);
  hPt->Scale(lumi/Nsim1);
  hphih->Scale(lumi/Nsim1);

  c00->cd(1);
  hx->DrawClone();
  c00->cd(2);
  hQ2->DrawClone();
  c00->cd(3);
  hz->DrawClone();
  c00->cd(4);
  hPt->DrawClone();
  c00->cd(5);
  hphih->DrawClone();

  hx->Reset();
  hQ2->Reset();
  hz->Reset();
  hPt->Reset();
  hphih->Reset();
  hx->SetLineColor(2);
  hQ2->SetLineColor(2);
  hz->SetLineColor(2);
  hPt->SetLineColor(2);
  hphih->SetLineColor(2);

  for (Long64_t i = 0; i < Nsim2; i++){
    if (i%100 == 0) cout << i << endl;
    weight = mysidis.GibbsSampler(0, 1);
    acc = 0;
    lp = mysidis.GetLorentzVector("lp");
    Ph = mysidis.GetLorentzVector("Ph");
    x = mysidis.GetVariable("x");
    Q2 = mysidis.GetVariable("Q2");
    z = mysidis.GetVariable("z");
    Pt = mysidis.GetVariable("Pt");
    W = mysidis.GetVariable("W");
    Wp = mysidis.GetVariable("Wp");
    phih = mysidis.GetVariable("phih");

    if (Q2 < 1.0) continue;
    if (W < 2.3) continue;
    if (Wp < 1.6) continue;
    
    if (ic == 1)
      acc = GetAcceptance(lp, "e-") * GetAcceptance(Ph, "h+");
    else if (ic == -1)
      acc = GetAcceptance(lp, "e-") * GetAcceptance(Ph, "h-");
    
    if (acc > 0){
      hx->Fill(x, weight * acc);
      hQ2->Fill(Q2, weight * acc);
      hz->Fill(z, weight * acc);
      hPt->Fill(Pt, weight * acc);
      hphih->Fill(phih, weight * acc);
    }
  }

  hx->Scale(lumi/Nsim2);
  hQ2->Scale(lumi/Nsim2);
  hz->Scale(lumi/Nsim2);
  hPt->Scale(lumi/Nsim2);
  hphih->Scale(lumi/Nsim2);

  c00->cd(1);
  hx->DrawClone("same");
  c00->cd(2);
  hQ2->DrawClone("same");
  c00->cd(3);
  hz->DrawClone("same");
  c00->cd(4);
  hPt->DrawClone("same");
  c00->cd(5);
  hphih->DrawClone("same");

  c00->Print("compare_gibbs.pdf");
      
  return 0;
}
