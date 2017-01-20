#include "Lsidis.h"

using namespace std;

int main(int argc, char * argv[]){
  gRandom->SetSeed(0);
  
  TLorentzVector P(0, 0, 0, 0.938272);//Nucleon
  TLorentzVector l(0, 0, 11.0, 11.0);//electron
  Lsidis mysidis;
  mysidis.SetNucleus(1, 0);//proton number, neutron number
  mysidis.SetInitialState(l, P);//set initial state
  mysidis.SetHadron("pi+");//set detected hadron: pi+, pi-, pi0, K+, K-, K0, p
  mysidis.SetPDFset("CJ15lo");//choose a PDF set
  mysidis.SetFFset("DSSFFlo");//choose an FF set
  
  TLorentzVector lp, Ph;//scattered electron and detected hadron
  double weight;//weight of event
  
  double Xmin[6] = {0.05, 1.0, 0.3, 0.0, -M_PI, -M_PI};//generator kinematic range x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.75, 8.0, 0.7, 1.0, M_PI, M_PI};
  
  mysidis.SetRange(Xmin, Xmax);
  
  for (int i = 0; i < 100; i++){
    weight = mysidis.GenerateEvent(0, 1);//generate an event
    if (weight > 0){
      lp = mysidis.GetLorentzVector("lp");
      Ph = mysidis.GetLorentzVector("Ph");
      cout << i << endl;
    }
  }
  
  return 0;
}
