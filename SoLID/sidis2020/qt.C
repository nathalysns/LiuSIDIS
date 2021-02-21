#include "Lsidis.h"
#include "SoLID_SIDIS_3He.h"

using namespace std;

int main(const int argc, const char * argv[]){

  Long64_t Nsim = 10000000;

  Lsidis sidis;

  TLorentzVector l(0, 0, 11.0, 11.0);
  TLorentzVector P(0, 0, 0, 0.938272);

  sidis.SetNucleus(2.0, 1.0);
  sidis.SetHadron("pi+");
  sidis.ChangeTMDpars(0.57, 0.12);

  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("DSSFFlo");
  double lumi = 1.0e+10 * pow(0.197327, 2);
  double time1 = 48.0 * 24.0 * 3600;
  double time2 = 21.0 * 24.0 * 3600;
  double Xmin[6] = {0.04, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.7, 10.0, 0.7, 2.0, M_PI, M_PI};

  sidis.SetRange(Xmin, Xmax);

  double weight = 0.0;
  double acc = 0.0;
  double x, Q2, qt;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);

  TFile * fs = new TFile("phasespace/qt-3he.root", "RECREATE");
  TH2D * h0 = new TH2D("11GeV_qt<Q", "", 100, 0.0, 1.0, 100, 0.0, 10.0);
  TH2D * h1 = new TH2D("11GeV_qt<0.33Q", "", 100, 0.0, 1.0, 100, 0.0, 10.0);
  TH2D * h2 = new TH2D("11GeV_qt<0.2Q", "", 100, 0.0, 1.0, 100, 0.0, 10.0);
  TH2D * g0 = new TH2D("8.8GeV_qt<Q", "", 100, 0.0, 1.0, 100, 0.0, 10.0);
  TH2D * g1 = new TH2D("8.8GeV_qt<0.33Q", "", 100, 0.0, 1.0, 100, 0.0, 10.0);
  TH2D * g2 = new TH2D("8.8GeV_qt<0.2Q", "", 100, 0.0, 1.0, 100, 0.0, 10.0);

  h0->SetDirectory(fs);
  h1->SetDirectory(fs);
  h2->SetDirectory(fs);
  g0->SetDirectory(fs);
  g1->SetDirectory(fs);
  g2->SetDirectory(fs); 
  
  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 1);
    if (i%1000000 == 0) cout << i / 100000 << "%" << endl;
    if (weight > 0) {
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp) * GetAcceptance_hadron(Ph, "pi+");
      if (acc > 0){
	x = sidis.GetVariable("x");
	Q2 = sidis.GetVariable("Q2");
	qt = sidis.GetVariable("Pt") / sidis.GetVariable("z");
	h0->Fill(x, Q2, weight * acc);
	if (qt < 0.33 * sqrt(Q2))
	  h1->Fill(x, Q2, weight * acc);
	if (qt < 0.2 * sqrt(Q2))
	  h2->Fill(x, Q2, weight * acc);
      }
    }
  }

  h0->Scale(lumi*time1/Nsim);
  h1->Scale(lumi*time1/Nsim);
  h2->Scale(lumi*time1/Nsim);

  l.SetXYZT(0, 0, 8.8, 8.8);
  sidis.SetInitialState(l, P);

  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 1);
    if (i%1000000 == 0) cout << i / 100000 << "%" << endl;
    if (weight > 0) {
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = GetAcceptance_e(lp) * GetAcceptance_hadron(Ph, "pi+");
      if (acc > 0){
	x = sidis.GetVariable("x");
	Q2 = sidis.GetVariable("Q2");
	qt = sidis.GetVariable("Pt") / sidis.GetVariable("z");
	g0->Fill(x, Q2, weight * acc);
	if (qt < 0.33 * sqrt(Q2))
	  g1->Fill(x, Q2, weight * acc);
	if (qt < 0.2 * sqrt(Q2))
	  g2->Fill(x, Q2, weight * acc);
      }
    }
  }

  g0->Scale(lumi*time2/Nsim);
  g1->Scale(lumi*time2/Nsim);
  g2->Scale(lumi*time2/Nsim);

  fs->Write();

  return 0;
}
  
  
