#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

using namespace std;

const LHAPDF::PDF * xf = LHAPDF::mkPDF("CJ15lo", 0);
const LHAPDF::PDF * zDpion = LHAPDF::mkPDF("DSSFFlo", 211);
const LHAPDF::PDF * zDkaon = LHAPDF::mkPDF("DSSFFlo", 321);

const double Mp = 0.938272;
const double kt2 = 0.57;
const double pt2 = 0.12;

double Nuv, auv, buv;
double Ndv, adv, bdv;
double Nub, Ndb, M1;

double xf1(const int flavor, const double x, const double Q2, const char * target){
  double Np = 1.0;
  double Nn = 0.0;
  if (strcmp(target, "neutron") == 0){
    Np = 0.0;
    Nn = 1.0;
  }
  else if (strcmp(target, "deuteron") == 0){
    Np = 0.5;
    Nn = 0.5;
  }
  else if (strcmp(target, "helium3") == 0){
    Np = 0.67;
    Nn = 0.33;
  }
  if (abs(flavor) == 1 || abs(flavor) == 2){
    return Np * xf->xfxQ2(flavor, x, Q2) + Nn * xf->xfxQ2(2/flavor, x, Q2);
  }
  else {
    return (Np + Nn) * xf->xfxQ2(flavor, x, Q2);
  }
}

double zD1(const int flavor, const double z, const double Q2, const char * hadron){
  if (strcmp(hadron, "pi+") == 0 || strcmp(hadron, "h+") == 0)
    return zDpion->xfxQ2(flavor, z, Q2);
  else if (strcmp(hadron, "pi-") == 0 || strcmp(hadron, "h-") == 0)
    return zDpion->xfxQ2(-flavor, z, Q2);
  else if (strcmp(hadron, "k+") == 0)
    return zDkaon->xfxQ2(flavor, z, Q2);
  else if (strcmp(hadron, "k-") == 0)
    return zDkaon->xfxQ2(-flavor, z, Q2);
  else
    return 0;
}

double FUUT(const double * var, const char * target, const char * hadron){
  //var: Q2, x, y, z, pT
  double Q2 = var[0];
  double x = var[1];
  double z = var[3];
  double pT = var[4];
  double Pt2 = pt2 + z * z * kt2;
  double eu2 = pow(2.0 / 3.0, 2);
  double ed2 = pow(-1.0 / 3.0, 2);
  double coll = eu2 * (xf1(2, x, Q2, target) * zD1(2, z, Q2, hadron) + xf1(-2, x, Q2, target) * zD1(-2, z, Q2, hadron))
    + ed2 * (xf1(1, x, Q2, target) * zD1(1, z, Q2, hadron) + xf1(-1, x, Q2, target) * zD1(-1, z, Q2, target) + xf1(3, x, Q2, target) * zD1(3, z, Q2, hadron)+ xf1(-3, x, Q2, target) * zD1(-3, z, Q2, hadron));
  double result = coll / z * exp(-pT * pT / Pt2) / (M_PI * Pt2);
  return result;
}
  
double Nxfx(const int flavor, const double x, const double Q2, const char * target){
  if (abs(flavor) > 2) return 0;
  double uv = Nuv * pow(x, auv) * pow(1.0 - x, buv) * pow(auv + buv, auv + buv) / pow(auv, auv) / pow(buv, buv) * (xf->xfxQ2(2, x, Q2) - xf->xfxQ2(-2, x, Q2)) / x;
  double dv = Ndv * pow(x, adv) * pow(1.0 - x, bdv) * pow(adv + bdv, adv + bdv) / pow(adv, adv) / pow(bdv, bdv) * (xf->xfxQ2(1, x, Q2) - xf->xfxQ2(-1, x, Q2)) / x;
  double ub = Nub * xf->xfxQ2(-2, x, Q2) / x;
  double db = Ndb * xf->xfxQ2(-1, x, Q2) / x;
  if (strcmp(target, "proton") == 0){
    if (flavor == 2) return uv + ub;
    else if (flavor == 1) return dv + db;
    else if (flavor == -2) return ub;
    else if (flavor == -1) return db;
  }
  else if (strcmp(target, "neutron") == 0){
    if (flavor == 2) return dv + db;
    else if (flavor == 1) return uv + ub;
    else if (flavor == -2) return db;
    else if (flavor == -1) return ub;
  }
  else if (strcmp(target, "deuteron") == 0){
    if (flavor == 2 || flavor == 1)
      return 0.5 * (uv + ub + dv + db);
    else if (flavor == -2 || flavor == -1)
      return 0.5 * (ub + db);
  }
  else if (strcmp(target, "helium3") == 0){
    if (flavor == 2) return 0.86 * (dv + db) / 3;
    else if (flavor == 1) return 0.86 * (uv + ub) / 3;
    else if (flavor == -2) return 0.86 * db / 3;
    else if (flavor == -1) return 0.86 * ub / 3;
  }
  return 0;
}

double f1T(const int flavor, const double x, const double kt, const double Q2){
  double ks2 = M1 * M1 * kt2 / (M1 * M1 + kt2);
  return -Mp / M1 * sqrt(2.0 * M_E) * Nxfx(flavor, x, Q2, "proton") * exp(-kt * kt / ks2) / (M_PI * kt2);
}

double f1T1(const int flavor, const double x, const double Q2){
  double ks2 = M1 * M1 * kt2 / (M1 * M1 + kt2);
  return -sqrt(2.0 * M_E) / (2.0 * Mp * M1) * Nxfx(flavor, x, Q2, "proton") * pow(ks2, 2) / kt2;
}

double AUTSivers(const double * var, const char * target, const char * hadron){
  //var: Q2, x, y, z, pT
  double Q2 = var[0];
  double x = var[1];
  double z = var[3];
  double pT = var[4];
  double ks2 = M1 * M1 * kt2 / (M1 * M1 + kt2);
  double prefactor = (z * z * kt2 + pt2) * pow(ks2, 2) / pow(z * z * ks2 + pt2, 2) / kt2 * exp(-pT * pT * (z * z * (kt2 - ks2)) / (z * z * ks2 + pt2) / (z * z * kt2 + pt2)) * sqrt(2.0 * M_E) * z * pT / M1;
  double eu2 = pow(2.0 / 3.0, 2);
  double ed2 = pow(-1.0/ 3.0, 2);
  double Numerator = eu2 * (Nxfx(2, x, Q2, target) * zD1(2, z, Q2, hadron) + Nxfx(-2, x, Q2, target) * zD1(-2, z, Q2, hadron)) + ed2 * (Nxfx(1, x, Q2, target) * zD1(1, z, Q2, hadron) + Nxfx(-1, x, Q2, target) * zD1(-1, z, Q2, hadron));
  double Denominator = eu2 * (xf1(2, x, Q2, target) * zD1(2, z, Q2, hadron) + xf1(-2, x, Q2, target) * zD1(-2, z, Q2, hadron))
    + ed2 * (xf1(1, x, Q2, target) * zD1(1, z, Q2, hadron) + xf1(-1, x, Q2, target) * zD1(-1, z, Q2, target) + xf1(3, x, Q2, target) * zD1(3, z, Q2, hadron)+ xf1(-3, x, Q2, target) * zD1(-3, z, Q2, hadron));
  double result = prefactor * Numerator / Denominator * x;
  return result;
}
 

#endif
