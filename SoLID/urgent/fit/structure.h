#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

using namespace std;

const LHAPDF::PDF * xf = LHAPDF::mkPDF("CJ15lo", 0);
const LHAPDF::PDF * zDp = LHAPDF::mkPDF("DSSFFlo", 211);
const LHAPDF::PDF * zDm = LHAPDF::mkPDF("DSSFFlo", 1211);
const LHAPDF::PDF * zD0 = LHAPDF::mkPDF("DSSFFlo", 111);

double xPDF(const int flavor, const double x, const double Q2, const double Np, const double Nn){
  if (abs(flavor) == 1 || abs(flavor) == 2)
    return Np * xf->xfxQ2(flavor, x, Q2) + Nn * xf->xfxQ2(2 / flavor, x, Q2);
  else
    return (Np + Nn) * xf->xfxQ2(flavor, x, Q2);
}

double f1(const int flavor, const double x, const double Q2, const char * target){
  if (strcmp(target, "proton") == 0)
    return xPDF(flavor, x, Q2, 1.0, 0.0) / x;
  if (strcmp(target, "neutron") == 0)
    return xPDF(flavor, x, Q2, 0.0, 1.0) / x;
  if (strcmp(target, "deuteron") == 0)
    return xPDF(flavor, x, Q2, 0.5, 0.5) / x;
  return 0;
}

double D1(const int flavor, const double z, const double Q2, const char * hadron){
  if (strcmp(hadron, "pi+") == 0)
    return zDp->xfxQ2(flavor, z, Q2) / z;
  if (strcmp(hadron, "pi-") == 0)
    return zDm->xfxQ2(flavor, z, Q2) / z;
  if (strcmp(hadron, "pi0") == 0)
    return zD0->xfxQ2(flavor, z, Q2) / z;
  return 0;
}

double FUUT(const double * var, const char * target, const char * hadron){
  //var: Q2, x, y, z, pT
  double Q2 = var[0];
  double x = var[1];
  double z = var[3];
  double pT = var[4];
  double kt2 = 0.604;
  double pt2 = 0.114;
  double PT2 = pt2 + z * z * kt2;
  double fu = f1(2, x, Q2, target);
  double fub = f1(-2, x, Q2, target);
  double fd = f1(1, x, Q2, target);
  double fdb = f1(-1, x, Q2, target);
  double fs = f1(3, x, Q2, target);
  double fsb = f1(-3, x, Q2, target);
  double fc = f1(4, x, Q2, target);
  double fcb = f1(-4, x, Q2, target);
  double fb = f1(5, x, Q2, target);
  double fbb = f1(-5, x, Q2, target);
  double Du = D1(2, z, Q2, hadron);
  double Dub = D1(-2, z, Q2, hadron);
  double Dd = D1(1, z, Q2, hadron);
  double Ddb = D1(-1, z, Q2, hadron);
  double Ds = D1(3, z, Q2, hadron);
  double Dsb = D1(-3, z, Q2, hadron);
  double Dc = D1(4, z, Q2, hadron);
  double Dcb = D1(-4, z, Q2, hadron);
  double Db = D1(5, z, Q2, hadron);
  double Dbb = D1(-5, z, Q2, hadron);
  double col = pow(2.0 / 3.0, 2) * (fu * Du + fub * Dub + fc * Dc + fcb * Dcb) + pow(1.0 / 3.0, 2) * (fd * Dd + fdb * Ddb + fs * Ds + fsb * Dsb + fb * Db + fbb * Dbb);
  double result = x * col * exp(- pT * pT / PT2) / (M_PI * PT2);
  return result;
}

double xTransversity(const int flavor, const double x, const double Q2, const double Np, const double Nn, const double * par){
  if (flavor == 1 || flavor == 2){
    double Nu = par[0];
    double au = par[1];
    double bu = par[2];
    double Nd = par[3];
    double ad = par[4];
    double bd = par[5];
    double xhu = Nu * pow(x, au) * pow(1.0 - x, bu) * xf->xfxQ2(2, x, Q2);
    double xhd = Nd * pow(x, ad) * pow(1.0 - x, bd) * xf->xfxQ2(1, x, Q2);
    if (flavor == 2)
      return Np * xhu + Nn * xhd;
    else
      return Np * xhd + Nn * xhu;
  }
  else
    return 0;
}

double h1(const int flavor, const double x, const double Q2, const char * target, const double * par){
  if (strcmp(target, "proton") == 0)
    return xTransversity(flavor, x, Q2, 1.0, 0.0, par) / x;
  if (strcmp(target, "neutron") == 0)
    return xTransversity(flavor, x, Q2, 0.0, 1.0, par) / x;
  if (strcmp(target, "deuteron") == 0)
    return xTransversity(flavor, x, Q2, 0.5, 0.5, par) / x;
  return 0;
}

double H1(const int flavor, const double z, const double Q2, const char * hadron){
  double Nfav = 0.49;
  double Ndis = -1.0;
  double c = 1.06;
  double d = 0.07;
  double Mh = sqrt(1.5);
  if ( (flavor == 2 && strcmp(hadron, "pi+") == 0) || (flavor == 1 && strcmp(hadron, "pi-") == 0) ){
    return sqrt(2.0 * M_E) * 0.14 * Mh / (Mh * Mh + 0.2) * Nfav * pow(z, c) * pow(1.0 - z, d) * pow(c + d, c + d) / pow(c, c) / pow(d, d) * zDp->xfxQ2(2, z, Q2);
  }
  else if ( (flavor == 1 && strcmp(hadron, "pi+") == 0) || (flavor == 2 && strcmp(hadron, "pi+") == 0) ){
    return sqrt(2.0 * M_E) * 0.14 * Mh / (Mh * Mh + 0.2) * Ndis * pow(z, c) * pow(1.0 - z, d) * pow(c + d, c + d) / pow(c, c) / pow(d, d) * zDp->xfxQ2(1, z, Q2);
  }
  else
    return 0;
}

double FUTCollins(const double * var, const char * target, const char * hadron, const double * par){
  //var: Q2, x, y, z, pT
  double Q2 = var[0];
  double x = var[1];
  double z = var[3];
  double pT = var[4];
  double kt2 = par[6];
  double pt2 = 1.5 * 0.2 / (1.5 + 0.2);
  double PT2 = pt2 + z * z * kt2;
  double factor = pt2 * pT / (0.14 * z * PT2);
  double hu = h1(2, x, Q2, target, par);
  double hd = h1(1, x, Q2, target, par);
  double Hu = H1(2, z, Q2, hadron);
  double Hd = H1(1, z, Q2, hadron);
  double col = pow(2.0 / 3.0, 2) * hu * Hu + pow(1.0 / 3.0, 2) * hd * Hd;
  double result = x * factor * col * exp(- pT * pT / PT2) / (M_PI * PT2);
  return result;
}

double AUTCollins(const double * var, const char * target, const char * hadron, const double * par){
  //var: Q2, x, y, z, pT
  double Q = sqrt(var[0]);
  double x = var[1];
  double y = var[2];
  double g = 2.0 * x * 0.939 / Q;
  double epsilon = (1.0 - y - 0.25 * g * g * y * y) / (1.0 - y + 0.5 * y * y + 0.25 * g * g * y * y);
  double result = epsilon * FUTCollins(var, target, hadron, par) / FUUT(var, target, hadron);
  return result;
}
 


#endif
