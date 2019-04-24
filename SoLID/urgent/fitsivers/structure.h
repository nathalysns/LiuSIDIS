#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

#include "Math/WrappedParamFunction.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/GSLIntegrator.h"

using namespace std;

const LHAPDF::PDF * xf = LHAPDF::mkPDF("CJ15lo", 0);
const LHAPDF::PDF * zDp = LHAPDF::mkPDF("DSSFFlo", 211);


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
  if (strcmp(hadron, "pi+") == 0 || strcmp(hadron, "h+") == 0)
    return zDp->xfxQ2(flavor, z, Q2) / z;
  if (strcmp(hadron, "pi-") == 0 || strcmp(hadron, "h-") == 0)
    return zDp->xfxQ2(-flavor, z, Q2) / z;
  return 0;
}

double FUUT(const double * var, const char * target, const char * hadron){
  //var: Q2, x, y, z, pT
  double Q2 = var[0];
  double x = var[1];
  double z = var[3];
  double pT = var[4];
  double kt2 = 0.25;
  double pt2 = 0.20;
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

double xf1T1(const int flavor, const double x, const double Q2, const double Np, const double Nn, const double * par){
  if (abs(flavor) == 1 || abs(flavor) == 2){
    double Nu = par[0];
    double au = par[1];
    double bu = par[2];
    double cu = par[3];
    double Nd = par[4];
    double ad = par[5];
    double bd = par[6];
    double cd = par[7];
    double Nub = par[8];
    double Ndb = par[9];
    double xf1T1u = Nu * (1.0 + cu * x) * pow(x, au) * pow(1.0 - x, bu) * pow(au + bu, au + bu) / pow(au, au) / pow(bu, bu) * xf->xfxQ2(2, x, Q2);
    double xf1T1d = Nd * (1.0 + cd * x) * pow(x, ad) * pow(1.0 - x, bd) * pow(ad + bd, ad + bd) / pow(ad, ad) / pow(bd, bd) * xf->xfxQ2(1, x, Q2);
    double xf1T1ub = Nub * xf->xfxQ2(-2, x, Q2);
    double xf1T1db = Ndb * xf->xfxQ2(-1, x, Q2);
    if (flavor == 2)
      return Np * xf1T1u + Nn * xf1T1d;
    else if (flavor == 1)
      return Np * xf1T1d + Nn * xf1T1u;
    else if (flavor == -2)
      return Np * xf1T1ub + Nn * xf1T1db;
    else
      return Np * xf1T1db + Nn * xf1T1ub;
  }
  else
    return 0;
}

double f1T1(const int flavor, const double x, const double Q2, const char * target, const double * par){
  if (strcmp(target, "proton") == 0)
    return xf1T1(flavor, x, Q2, 1.0, 0.0, par) / x;
  if (strcmp(target, "neutron") == 0)
    return xf1T1(flavor, x, Q2, 0.0, 1.0, par) / x;
  if (strcmp(target, "deuteron") == 0)
    return xf1T1(flavor, x, Q2, 0.5, 0.5, par) / x;
  return 0;
}

double FUTSivers(const double * var, const char * target, const char * hadron, const double * par){
  //var: Q2, x, y, z, pT
  double Q2 = var[0];
  double x = var[1];
  double z = var[3];
  double pT = var[4];
  double kt2 = par[10];
  double pt2 = 0.20;
  double PT2 = pt2 + z * z * kt2;
  double factor = -2.0 * z * 0.939 * pT / PT2;
  double f1T1u = f1T1(2, x, Q2, target, par);
  double f1T1d = f1T1(1, x, Q2, target, par);
  double f1T1ub = f1T1(-2, x, Q2, target, par);
  double f1T1db = f1T1(-1, x, Q2, target, par);
  double D1u = D1(2, z, Q2, hadron);
  double D1d = D1(1, z, Q2, hadron);
  double D1ub = D1(-2, z, Q2, hadron);
  double D1db = D1(-1, z, Q2, hadron);
  double col = pow(2.0 / 3.0, 2) * (f1T1u * D1u + f1T1ub * D1ub) + pow(1.0 / 3.0, 2) * (f1T1d * D1d + f1T1db * D1db);
  double result = x * factor * col * exp(- pT * pT / PT2) / (M_PI * PT2);
  return result;
}

double AUTSivers(const double * var, const char * target, const char * hadron, const double * par){
  //var: Q2, x, y, z, pT
  double result = FUTSivers(var, target, hadron, par) / FUUT(var, target, hadron);
  return result;
}


#endif
