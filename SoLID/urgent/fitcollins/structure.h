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

const LHAPDF::PDF * xf = LHAPDF::mkPDF("CT14lo", 0);
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
  if (strcmp(hadron, "pi+") == 0)
    return zDp->xfxQ2(flavor, z, Q2) / z;
  if (strcmp(hadron, "pi-") == 0)
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

double xTransversity(const int flavor, const double x, const double Q2, const double Np, const double Nn, const double * par){
  if (flavor == 1 || flavor == 2){
    double Nu = par[0];
    double Nd = par[1];
    double a = par[2];
    double b = par[3];
    double c = par[4];
    double xhu = Nu * (1.0 + c * sqrt(x)) * pow(x, a) * pow(1.0 - x, b) * pow(a + b, a + b) / pow(a, a) / pow(b, b) * (xf->xfxQ2(2, x, Q2) - xf->xfxQ2(-2, x, Q2));
    double xhd = Nd * (1.0 + c * sqrt(x)) * pow(x, a) * pow(1.0 - x, b) * pow(a + b, a + b) / pow(a, a) / pow(b, b) * (xf->xfxQ2(1, x, Q2) - xf->xfxQ2(-1, x, Q2));
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
  double Nfav = 1.0;
  double Ndis = -1.0;
  double c = -2.36;
  double d = 2.12;
  double Mh = sqrt(0.67);
  double factor = sqrt(2.0 * M_E) * 0.14 * Mh / (Mh * Mh + 0.20);
  if ( (flavor == 2 && strcmp(hadron, "pi+") == 0) || (flavor == 1 && strcmp(hadron, "pi-") == 0) ){
    return factor * Nfav * ((1.0 - c - d) + c * z + d * z * z) * zDp->xfxQ2(2, z, Q2);
  }
  else if ( (flavor == 1 && strcmp(hadron, "pi+") == 0) || (flavor == 2 && strcmp(hadron, "pi-") == 0) ){
    return factor * Ndis * ((1.0 - c - d) + c * z + d * z * z) * zDp->xfxQ2(1, z, Q2);
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
  double kt2 = par[5];
  double pt2 = 0.67 * 0.20 / (0.67 + 0.20);
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

double h1u(double x, void * par){
  //double x = var[0];
  double Q2 = 2.4;
  return xTransversity(2, x, Q2, 1.0, 0.0, (double *) par) / x;
}

double h1d(double x, void * par){
  //double x = var[0];
  double Q2 = 2.4;
  return xTransversity(1, x, Q2, 1.0, 0.0, (double *) par) / x;
}

double gtu(double * par){
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-8, 1e-4);
  ig.SetFunction(&h1u, par);
  double result = ig.Integral(1e-5, 1.0);
  return result;
}

double gtd(double * par){
  ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, 1e-8, 1e-4);
  ig.SetFunction(&h1d, par);
  double result = ig.Integral(1e-5, 1.0);
  return result;
}
 


#endif
