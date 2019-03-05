#ifndef _EICCACCEPTANCE_H_
#define _EICCACCEPTANCE_H_

#include "TLorentzVector.h"

namespace EICC{

  double Rfactor0 = 1.0e5;
  
  double GetAcceptance_e(const TLorentzVector p){
    return 1;
  }
  
  double GetAcceptance_pip(const TLorentzVector p){
    return 1;
  }

  double GetAcceptance_pim(const TLorentzVector p){
    return 1;
  }

  double GetAcceptance_kp(const TLorentzVector p){
    return 1;
  }

  double GetAcceptance_km(const TLorentzVector p){
    return 1;
  }

  double GetAcceptance_hadron(const TLorentzVector p, const char * hadron){
    if (strcmp(hadron, "pi+") == 0) return GetAcceptance_pip(p);
    else if (strcmp(hadron, "pi-") == 0) return GetAcceptance_pim(p);
    else if (strcmp(hadron, "K+") == 0) return GetAcceptance_kp(p);
    else if (strcmp(hadron, "K-") == 0) return GetAcceptance_km(p);
    else return 0;
  }

}

#endif
