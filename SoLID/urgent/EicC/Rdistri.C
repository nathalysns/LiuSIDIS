#include "sidis-region.h"
#include "eiccacceptance.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"


using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./Rdistri <opt>" << endl;
    cout << "./Rdistri <1> <hadron>" << endl;
    return 0;
  }
  
  const int opt = atoi(argv[1]);

  gRandom->SetSeed(0);

  const double Mp = 0.938272;

  TLorentzVector l(0, 0, 3.5, 3.5);
  TLorentzVector P1(0, 0, -sqrt(pow(20.0, 2) - Mp*Mp), 20.0);
  TLorentzVector P2(0, 0, -sqrt(pow(20.0*2/3, 2) - Mp*Mp), 20.0*2/3);
  TLorentzVector lp, Ph;
  double weight, acc, R1;
  double Xmin[6] = {0.001, 1.0, 0.3, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {0.999, 50.0, 0.7, 2.0, M_PI, M_PI};
  Long64_t Nsim = 10000000;
  double lumi = 1e+7 * pow(0.197327, 2);
  double time = 100.0 * 24.0 * 3600.0;

  Lsidis sidis1;
  Lsidis sidis2;

  sidis1.SetNucleus(1, 0);//proton
  sidis1.SetInitialState(l, P1);
  sidis1.SetPDFset("CJ15lo");
  sidis1.ChangeTMDpars(0.57, 0.12);
  sidis1.SetRange(Xmin, Xmax);

  sidis2.SetNucleus(0.67, 0.33);//helium-3
  sidis2.SetInitialState(l, P2);
  sidis2.SetPDFset("CJ15lo");
  sidis2.ChangeTMDpars(0.57, 0.12);
  sidis2.SetRange(Xmin, Xmax);
  
  if (opt == 0){//test
    return 0;
  }

  if (opt == 1){//R1 analysis
    sidis1.SetHadron(argv[2]);
    sidis1.SetFFset("DSSFFlo");
    sidis2.SetHadron(argv[2]);
    sidis2.SetFFset("DSSFFlo");

    TFile * fs = new TFile("R1.root", "RECREATE");
    TH1D * h1 = new TH1D("h1", "", 200, 0.0, 10.0);
    h1->SetDirectory(fs);
    TH1D * h2 = new TH1D("h2", "", 200, 0.0, 10.0);
    h2->SetDirectory(fs);
    TH2D * zPt1 = new TH2D("zPt1", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1->SetDirectory(fs);
    TH2D * zPt1cut = new TH2D("zPt1cut", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1cut->SetDirectory(fs);
    TH2D * zPt2 = new TH2D("zPt2", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1->SetDirectory(fs);
    TH2D * zPt2cut = new TH2D("zPt2cut", "", 200, 0.0, 1.0, 200, 0.0, 2.0);
    zPt1cut->SetDirectory(fs);

    int proc = 0;
    for (Long64_t i = 0; i < Nsim; i++){//proton
      if (i * 100 / Nsim > proc){
	proc = i * 100 / Nsim;
	cout << proc << "%" << endl;
      }
      weight = sidis1.GenerateEvent(0, 1);
      if (weight > 0){
	if (sidis1.GetVariable("W") < 2.3) continue;
	if (sidis1.GetVariable("Wp") < 1.6) continue;
	if (sidis1.GetVariable("Q2") < 1.0) continue;
	if (sidis1.GetVariable("z") < 0.3 || sidis1.GetVariable("z") > 0.7) continue;
	lp = sidis1.GetLorentzVector("lp");
	Ph = sidis1.GetLorentzVector("Ph");
	acc = EICC::GetAcceptance_e(lp) * EICC::GetAcceptance_hadron(Ph, argv[2]);
	if (acc > 0){
	  R1 = Calculate_R1(&sidis1);
	  h1->Fill(R1, weight * acc);
	  zPt1->Fill(sidis1.GetVariable("z"), sidis1.GetVariable("Pt"), weight * acc);
	  if (R1 < 0.4){
	    zPt1cut->Fill(sidis1.GetVariable("z"), sidis1.GetVariable("Pt"), weight * acc);
	  }
	}
      }

      weight = sidis2.GenerateEvent(0, 1);
      if (weight > 0){
	if (sidis2.GetVariable("W") < 2.3) continue;
	if (sidis2.GetVariable("Wp") < 1.6) continue;
	if (sidis2.GetVariable("Q2") < 1.0) continue;
	if (sidis2.GetVariable("z") < 0.3 || sidis2.GetVariable("z") > 0.7) continue;
	lp = sidis2.GetLorentzVector("lp");
	Ph = sidis2.GetLorentzVector("Ph");
	acc = EICC::GetAcceptance_e(lp) * EICC::GetAcceptance_hadron(Ph, argv[2]);
	if (acc > 0){
	  R1 = Calculate_R1(&sidis2);
	  h2->Fill(R1, weight * acc);
	  zPt2->Fill(sidis2.GetVariable("z"), sidis2.GetVariable("Pt"), weight * acc);
	  if (R1 < 0.4){
	    zPt2cut->Fill(sidis2.GetVariable("z"), sidis2.GetVariable("Pt"), weight * acc);
	  }
	}
      }
    }
    h1->Scale(lumi*time/Nsim);
    h2->Scale(lumi*time/Nsim);
    zPt1->Scale(lumi*time/Nsim);
    zPt1cut->Scale(lumi*time/Nsim);
    zPt2->Scale(lumi*time/Nsim);
    zPt2cut->Scale(lumi*time/Nsim);

    fs->Write();
    return 0;
  }
  
  

  return 0;
}
