#include "SoLID_SIDIS_3He.h"

using namespace std;

int main(){

  gRandom->SetSeed(0);
  
  Lsidis sidis_3he;
  Lsidis sidis_neutron;
  
  TLorentzVector l(0, 0, 8.8, 8.8);
  TLorentzVector P(0, 0, 0, 0.938272);

  double lumi = 1.0e+10 * pow(0.197327, 2);
  double Nsim = 1.0e6;

  double Xmin[6] = {0.03, 1.0, 0.3, 0.0, -M_PI, -M_PI};//x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.7, 10.0, 0.7, 1.8, M_PI, M_PI};//x, Q2, z, Pt, phih, phiS

  sidis_3he.SetNucleus(2, 1);
  sidis_3he.SetHadron("pi+");
  sidis_3he.SetInitialState(l, P);
  sidis_3he.SetPDFset("CT14lo");
  sidis_3he.SetFFset("DSSFFlo");
  sidis_3he.SetRange(Xmin, Xmax);

  sidis_neutron.SetNucleus(0, 1);
  sidis_neutron.SetHadron("pi+");
  sidis_neutron.SetInitialState(l, P);
  sidis_neutron.SetPDFset("CT14lo");
  sidis_neutron.SetFFset("DSSFFlo");
  sidis_neutron.SetRange(Xmin, Xmax);

  double weight = 0.0;
  double weight2 = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);

  TH1D * hx = new TH1D("hx", "", 10, 0.0, 1.0);
  TH1D * nx = new TH1D("nx", "", 10, 0.0, 1.0);

  double x;
  double acc;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i%10000000 == 0) cout << i << endl;
    weight = sidis_3he.GenerateEvent(0, 1);
    if (weight > 0){
      lp = sidis_3he.GetLorentzVector("lp");
      Ph = sidis_3he.GetLorentzVector("Ph");
      
      acc = GetAcceptance_e(lp) * GetAcceptance_pi(Ph);
      
      if (acc > 0){
	x = sidis_3he.GetVariable("x");
	sidis_neutron.SetFinalState(lp, Ph);
	sidis_neutron.CalculateVariables();
	//cout << x << " " << sidis_neutron.GetVariable("x") << endl;
	weight2 = sidis_neutron.GetEventWeight(0, 1);
	
	hx->Fill(x, weight * acc);
	nx->Fill(x, weight2 * acc);
      }
    }
  }

  hx->Scale(lumi/Nsim);
  nx->Scale(lumi/Nsim);

  cout << hx->Integral(1, -1) << endl;
  cout << nx->Integral(1, -1) << endl;

  return 0;
}
	

	



      

  
  

