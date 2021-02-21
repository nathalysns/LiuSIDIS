#include "Lsidis.h"
#include "SoLID_SIDIS_3He.h"

using namespace std;

double EICAcceptance(const TLorentzVector P, const char * part = "e"){
  double rap = log((P.E() + P.Pz()) / (P.E() - P.Pz())) / 2.0;
  if (rap > 3.5 || rap < -3.5) return 0;
  double p = P.P();
  if (strcmp(part, "e") == 0 || strcmp(part, "e-") == 0){
    if (rap < -1.0 && p > 18.0) return 0;
    if (rap > -1.0 && rap < 1.0 && p > 8.0) return 0;
    if (rap > 1.0 && p > 20.0) return 0;
  }
  if (strcmp(part, "pi+") == 0 || strcmp(part, "pi-") == 0){
    if (rap < -1.0 && p > 7.0) return 0;
    if (rap > -1.0 && rap < 1.0 && p > 9.0) return 0;
    if (rap > 1.0 && p > 50.0) return 0;
  }
  return 1.0;
}

int main(const int argc, const char * argv[]){

  Long64_t Nsim = 10000000;

  Lsidis sidis;

  TLorentzVector l(0, 0, -5.0, 5.0);
  TLorentzVector P(0, 0, sqrt(pow(41.0, 2) - pow(0.938272,2)), 41.0);

  sidis.SetNucleus(1.0, 0.0);
  sidis.SetHadron("pi+");
  sidis.ChangeTMDpars(0.57, 0.12);

  sidis.SetInitialState(l, P);
  sidis.SetPDFset("CJ15lo");
  sidis.SetFFset("JAM19FF_pion_nlo");
  double lumi = 10.0 * 1.0e+13 * pow(0.197327, 2);
  
  double Xmin[6] = {1.0e-4, 0.01, 0.2, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {1.0, 0.95, 0.8, 2.0, M_PI, M_PI};

  sidis.SetRange(Xmin, Xmax);

  double weight = 0.0;
  double acc = 0.0;
  double x, Q2, qt;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);

  double X[101];
  double Y[101];
  for (int i = 0; i < 101; i++){
    X[i] = pow(10.0, -4.0 + 4.0 / 100 * i);
    Y[i] = pow(10.0, -1.0 + 4.0 / 100 * i);
  }

  
  TFile * fs = new TFile("phasespace/qt-eic.root", "RECREATE");
  TH2D * h0 = new TH2D("qt<Q", "", 100, X, 100, Y);
  TH2D * h1 = new TH2D("qt<0.33Q", "", 100, X, 100, Y);
  TH2D * h2 = new TH2D("qt<0.2Q", "", 100, X, 100, Y);

  h0->SetDirectory(fs);
  h1->SetDirectory(fs);
  h2->SetDirectory(fs);
  
  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 0);
    if (i%1000000 == 0) cout << i / 100000 << "%" << endl;
    if (weight > 0) {
      if (sidis.GetVariable("Q2") < 1.0) continue;
      if (sidis.GetVariable("W") < sqrt(10.0)) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;
      lp = sidis.GetLorentzVector("lp");
      Ph = sidis.GetLorentzVector("Ph");
      acc = EICAcceptance(lp, "e") * EICAcceptance(Ph, "pi+");
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

  h0->Scale(lumi/Nsim);
  h1->Scale(lumi/Nsim);
  h2->Scale(lumi/Nsim);

  fs->Write();

  return 0;
}
  
  
