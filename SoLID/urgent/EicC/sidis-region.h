#ifndef _SIDIS_REGION_H_
#define _SIDIS_REGION_H_

#include "Lsidis.h"
#include "TF1.h"

double Maxwell(const double * v, const double * par = 0){
  return v[0] * v[0] * exp(-v[0] * v[0] / 0.02);
}

double xi_distri(const double * xi, const double * par = 0){
  return xi[0] * pow(1.0 - xi[0], 3);
}

double zeta_distri(const double * zeta, const double * par = 0){
  return zeta[0] * pow(1.0 - zeta[0], 2);
}

double kt_distri(const double * kt, const double * par = 0){
  return exp(- kt[0] / 0.4);
}

TF1 TF_Maxwell("Maxwell", Maxwell, 0.0, 1.0, 0);
TF1 TF_xi("xi", xi_distri, 0.0, 1.0, 0);
TF1 TF_zeta("zeta", zeta_distri, 0.0, 1.0, 0);
TF1 TF_kt("kt", kt_distri, 0.0, 1.5, 0);
  
double Calculate_zn(Lsidis * sidis){
  double Mp = 0.938272081;
  double Mh = sidis->GetVariable("Mh");
  double Pt = sidis->GetVariable("Pt");
  double xn = sidis->GetVariable("xn");
  double x = sidis->GetVariable("x");
  double Q2 = sidis->GetVariable("Q2");
  double z = sidis->GetVariable("z");
  double zn = xn * z / (2.0 * x) * (1.0 + sqrt(1.0 - 4.0 * Mp * Mp * (Mh * Mh + Pt * Pt) * x * x / (z * z * Q2 * Q2)));
  return zn;
}

double Calculate_R1(Lsidis * sidis){
  double Mh = sidis->GetVariable("Mh");
  double Pt = sidis->GetVariable("Pt");
  double xn = sidis->GetVariable("xn");
  double Q2 = sidis->GetVariable("Q2");
  double zn = Calculate_zn(sidis);
  double xi = TF_xi.GetRandom(xn, 1.0);
  double zeta = TF_zeta.GetRandom(zn, 1.0);
  double dkt = TF_kt.GetRandom();
  double ki = TF_Maxwell.GetRandom();
  double kf = TF_Maxwell.GetRandom();
  double angle = gRandom->Uniform(-M_PI, M_PI);
  double xnhat = xn / xi;
  double denominator = (zn * Q2) / (2.0 * xnhat) - pow(ki, 2) * xn * (Mh * Mh + Pt * Pt) / (2.0 * zn * Q2);
  double kft2 = pow(Pt / zeta, 2) + dkt * dkt - 2.0 * dkt * Pt / zeta * cos(angle);
  double numerator = (kft2 + kf * kf) * zeta / 2.0 + (Mh * Mh + Pt * Pt) / (2.0 * zeta) + Pt * Pt / zeta - dkt * Pt * cos(angle);
  return abs(numerator / denominator);
}
  
    




#endif
