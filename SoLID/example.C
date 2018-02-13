#include "Lsidis.h"

using namespace std;

int main(const int argc, const char * argv[]){

  Lsidis sidis;//Create an object

  TLorentzVector l(0, 0, 11.0, 11.0);//initial state lepton
  TLorentzVector P(0, 0, 0, 0.938272);//initial state nucleon
  sidis.SetInitialState(l, P);//set initial state kinematics
  
  sidis.SetNucleus(2, 1);//number of proton and number of neutron
  sidis.SetHadron("K+");//final state hadron: pi+, pi-, pi0, K+, K-, K0, p

  sidis.ChangeTMDpars(0.604, 0.131);//use this for kaon
  sidis.SetPDFset("CJ15lo");//we use this in our new fit
  sidis.SetFFset("DSSFFlo");//we use this in our new fit

  double lumi = 1.0e+10 * pow(0.197327, 2);//luminosity (in unit of GeV^2 s^-1), this value corresponds to 1e36 cm^-2 s^-1;
  //It is equivalent to use 3e36 cm^-2 s^-1 and sidis.SetNucleus(0.67, 0.33)
  double time = 3600.0;//in unit of s
  //only make sure the unit used for lumi * time cancels the unit for the cross section (I choose GeV^-2 for the cross section unit)

  double Xmin[6] = {0.0, 1.0, 0.3, 0.0, -M_PI, -M_PI};//lower limits of x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.7, 10.0, 0.7, 1.8, M_PI, M_PI};//upper limits of x, Q2, z, Pt, phih, phiS
  sidis.SetRange(Xmin, Xmax);//Set generator's kinematic range

  Long64_t Nsim = 100;//

  double weight = 0.0;//to save the weight
  TLorentzVector lp(0, 0, 0, 0);//to save scattered lepton
  TLorentzVector Ph(0, 0, 0, 0);//to save final state hadron

  for (Long64_t i = 0; i < Nsim; i++){
    weight = sidis.GenerateEvent(0, 1);//generate an event, it returns the diffrential cross section * generate volume, the unit is GeV^-2
    if (weight > 0){// nonzero cross section
      //if you want some kinematics cuts, e.g.
      if (sidis.GetVariable("W") < 2.3) continue;
      if (sidis.GetVariable("Wp") < 1.6) continue;

      lp = sidis.GetLorentzVector("lp");//Get the scattered lepton
      Ph = sidis.GetLorentzVector("Ph");//Get the final state hadron
      //Then you can compare them with Zhiwen's acceptance files
    }
  }
  
  return 0;
}
