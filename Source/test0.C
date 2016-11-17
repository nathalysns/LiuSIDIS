#include "../Header/Lsidis1.h"

using namespace std;

int main(){

  const double rad = M_PI / 180.0;

  TLorentzVector l(0, 0, 11.0, 11.0);
  TLorentzVector P(0, 0, 0, 0.938272);
  TLorentzVector lp, Ph;

  Lsidis mysidis;

  TLorentzVector boo(0.1, 0.15, -0.3, 1.0);
  //l.Boost(boo.BoostVector());
  //P.Boost(boo.BoostVector());

  mysidis.SetNucleus(2.0, 1.0);
  mysidis.SetHadron("pi+");
  mysidis.SetPDFset("CT14lo");
  mysidis.SetFFset("DSSFFlo");
  mysidis.SetInitialState(l, P);

  double Xmin[6] = {0.05, 0.05, 0.05, 0.0, -M_PI, -M_PI};
  double Xmax[6] = {1.0, 1.0, 1.0, 2.0, M_PI, M_PI};

  gRandom->SetSeed(0);
  mysidis.SetRange(Xmin, Xmax);
  for (int i = 0; i < 10; i++){
    mysidis.GenerateEvent(0, 0);
    lp = mysidis.GetLorentzVector("lp");
    Ph = mysidis.GetLorentzVector("Ph");
  }

  mysidis.SetVariables(0.192311, 0.665824, 0.359942, 0.290156, 2.00175, 0.334765);
  mysidis.CalculateFinalState();
  lp = mysidis.GetLorentzVector("lp");
  Ph = mysidis.GetLorentzVector("Ph");

  mysidis.SetVariables(0.132311, 0.665824, 0.359942, 0.290156, 2.00175, 0.334765);
  mysidis.CalculateFinalState();
  lp = mysidis.GetLorentzVector("lp");
  Ph = mysidis.GetLorentzVector("Ph");

  //cout << lp.E() << " " << lp.Theta() * 180.0 / M_PI << " " << lp.Phi() * 180 / M_PI << endl;
  cout << Ph.E() << " " << Ph.Theta() * 180.0 / M_PI << " " << Ph.Phi() * 180 / M_PI << endl;

  mysidis.SetFinalState(lp, Ph);
  mysidis.CalculateVariables();
  mysidis.Test();
  


  return 0;
}
