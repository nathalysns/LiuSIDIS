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
  double vz;
  
  Long64_t Nsim = 0;
  Long64_t Nrec = 0;
  double acc, accr;

  TH1D * hx = new TH1D("hx", "", 100, 0.0, 0.7);
  TH1D * hQ2 = new TH1D("hQ2", "", 100, 0.0, 10.0);
  TH1D * hz = new TH1D("hz", "", 100, 0.25, 0.75);
  TH1D * hPt = new TH1D("hPt", "", 100, 0.0, 1.6);
  TH1D * hphih = new TH1D("hphih", "", 100, -3.2, 3.2);

  TH1D * rx = new TH1D("rx", "", 100, 0.0, 0.7);
  TH1D * rQ2 = new TH1D("rQ2", "", 100, 0.0, 10.0);
  TH1D * rz = new TH1D("rz", "", 100, 0.25, 0.75);
  TH1D * rPt = new TH1D("rPt", "", 100, 0.0, 1.6);
  TH1D * rphih = new TH1D("rphih", "", 100, -3.2, 3.2);

  TH1D * rhx = new TH1D("rhx", "", 100, 0.0, 0.7);
  TH1D * rhQ2 = new TH1D("rhQ2", "", 100, 0.0, 10.0);
  TH1D * rhz = new TH1D("rhz", "", 100, 0.25, 0.75);
  TH1D * rhPt = new TH1D("rhPt", "", 100, 0.0, 1.6);
  TH1D * rhphih = new TH1D("rhphih", "", 100, -3.2, 3.2);

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
  hx->SetLineWidth(1.5);
  hQ2->SetLineColor(4);
  hQ2->SetLineWidth(1.5);
  hz->SetLineColor(4);
  hz->SetLineWidth(1.5);
  hPt->SetLineColor(4);
  hPt->SetLineWidth(1.5);
  hphih->SetLineColor(4);
  hphih->SetLineWidth(1.5);

  rx->SetFillColor(2);
  rQ2->SetFillColor(2);
  rz->SetFillColor(2);
  rPt->SetFillColor(2);
  rphih->SetFillColor(2);

  rx->SetLineColor(2);
  rx->SetLineWidth(0.05);
  rQ2->SetLineColor(2);
  rQ2->SetLineWidth(0.05);
  rz->SetLineColor(2);
  rz->SetLineWidth(0.05);
  rPt->SetLineColor(2);
  rPt->SetLineWidth(0.05);
  rphih->SetLineColor(2);
  rphih->SetLineWidth(0.05);


  rhx->GetXaxis()->SetTitle("x");
  rhx->SetLineColor(4);
  rhQ2->GetXaxis()->SetTitle("Q2");
  rhQ2->SetLineColor(4);
  rhz->GetXaxis()->SetTitle("z");
  rhz->SetLineColor(4);
  rhPt->GetXaxis()->SetTitle("Pt");
  rhPt->SetLineColor(4);
  rhphih->GetXaxis()->SetTitle("phih");
  rhphih->SetLineColor(4);
 
  double x, Q2, z, W, Wp, Pt, phih;

  while (Nrec < 100000){
    Nsim++;
    if (Nsim%1000000 == 0) cout << Nsim << endl;
    vz = gRandom->Uniform(-370, -330);
    weight = mysidis.GenerateEvent(0, 1);
    if (weight > 0){

      mysidis.CalculateRfactor();

      x = mysidis.GetVariable("Rfactor");
      if (true) cout << x << endl;

      acc = 0;
      accr = 0;
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
	Nrec++;
	accr = acc;
	if(lp.Angle(Ph.Vect()) < 2.5 / 180.0 * M_PI)
	  accr = acc / 10.0;
	
	hx->Fill(x, weight * acc);
	rx->Fill(x, weight * accr);

	hQ2->Fill(Q2, weight * acc);
	rQ2->Fill(Q2, weight * accr);
	
	hz->Fill(z, weight * acc);
	rz->Fill(z, weight * accr);
	
	hPt->Fill(Pt, weight * acc);
	rPt->Fill(Pt, weight * accr);

	hphih->Fill(phih, weight * acc);
	rphih->Fill(phih, weight * accr);
      }
    }
  }

  for (int i = 1; i < 370; i++){
    if (i <= hx->GetNbinsX() && hx->GetBinContent(i) > 0){
      rhx->SetBinContent(i, rx->GetBinContent(i)/hx->GetBinContent(i));
    }
    if (i <= hQ2->GetNbinsX() && hQ2->GetBinContent(i) > 0){
      rhQ2->SetBinContent(i, rQ2->GetBinContent(i)/hQ2->GetBinContent(i));
    }
    if (i <= hz->GetNbinsX() && hz->GetBinContent(i) > 0){
      rhz->SetBinContent(i, rz->GetBinContent(i)/hz->GetBinContent(i));
    }
    if (i <= hPt->GetNbinsX() && hPt->GetBinContent(i) > 0){
      rhPt->SetBinContent(i, rPt->GetBinContent(i)/hPt->GetBinContent(i));
    }
    if (i <= hphih->GetNbinsX() && hphih->GetBinContent(i) > 0){
      rhphih->SetBinContent(i, rphih->GetBinContent(i)/hphih->GetBinContent(i));
    }
  }
  
  rhx->GetXaxis()->SetTitle("x");
  rhx->GetXaxis()->SetTitleSize(0.05);
  rhx->GetXaxis()->SetLabelSize(0.05);
  rhx->GetYaxis()->SetLabelSize(0.05);
  rhx->GetYaxis()->SetRangeUser(0.6, 1.1);
  rhx->GetXaxis()->CenterTitle();
  rhQ2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  rhQ2->GetXaxis()->SetTitleSize(0.05);
  rhQ2->GetXaxis()->SetLabelSize(0.05);
  rhQ2->GetYaxis()->SetLabelSize(0.05);
  rhQ2->GetYaxis()->SetRangeUser(0.6, 1.1);
  rhQ2->GetXaxis()->CenterTitle();
  rhz->GetXaxis()->SetTitle("z");
  rhz->GetXaxis()->SetTitleSize(0.05);
  rhz->GetXaxis()->SetLabelSize(0.05);
  rhz->GetYaxis()->SetLabelSize(0.05);
  rhz->GetYaxis()->SetRangeUser(0.6, 1.1);
  rhz->GetXaxis()->CenterTitle();
  rhPt->GetXaxis()->SetTitle("P_{hT} (GeV)");
  rhPt->GetXaxis()->SetTitleSize(0.05);
  rhPt->GetXaxis()->SetLabelSize(0.05);
  rhPt->GetYaxis()->SetLabelSize(0.05);
  rhPt->GetYaxis()->SetRangeUser(0.6, 1.1);
  rhPt->GetXaxis()->CenterTitle();
  rhphih->GetXaxis()->SetTitle("#phi_{h} (rad)");
  rhphih->GetXaxis()->SetTitleSize(0.05);
  rhphih->GetXaxis()->SetLabelSize(0.05);
  rhphih->GetYaxis()->SetLabelSize(0.05);
  rhphih->GetYaxis()->SetRangeUser(0.6, 1.1);
  rhphih->GetXaxis()->CenterTitle();

  hx->Scale(lumi/Nsim);
  hQ2->Scale(lumi/Nsim);
  hz->Scale(lumi/Nsim);
  hPt->Scale(lumi/Nsim);
  hphih->Scale(lumi/Nsim);

  rx->Scale(lumi/Nsim);
  rQ2->Scale(lumi/Nsim);
  rz->Scale(lumi/Nsim);
  rPt->Scale(lumi/Nsim);
  rphih->Scale(lumi/Nsim);

  cout << hx->Integral(1, -1) << endl;
      
  return 0;
}
