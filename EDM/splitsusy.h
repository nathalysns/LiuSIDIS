#ifndef _SPLITSUSY_H_
#define _SPLITSUSY_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "gsl/gsl_sf_log.h"
#include "gsl/gsl_sf_dilog.h"

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/AlphaS.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TVector.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"

using namespace std;

const double MW = 80.385;
const double MZ = 91.1876;
const double VEV = 246.0;
const double MH = 125.09;
const double alpha = 1.0 / 137.0;
const double thetaW = asin(sqrt(0.2315));
const double Mu = 0.0022;
const double Md = 0.0047;
const double alphaS0 = 0.41 * M_PI;

const LHAPDF::PDF * pdf = LHAPDF::mkPDF("CJ15lo", 0);

TMatrixD X(2,2), U(2,2), V(2,2);
TVectorD MX(2);
TMatrixD Y(4,4), N(4,4);
TVectorD MY(4);
TMatrixD gDR(2,2), gDL(2,2), GR(2,2), GL(2,2);
TMatrixD CR(2,4), CL(2,4);

double dilog_re(const double re, const double im){//real part of the dilogrithm function
  double r = sqrt(re * re + im * im);
  double theta = acos(re / r);
  if (im < 0) theta = -theta;
  gsl_sf_result real, imaginary;
  gsl_sf_complex_dilog_e(r, theta, &real, &imaginary);
  return real.val;
}

double dilog_im(const double re, const double im){//imaginary part of the dilogrithm function
  double r = sqrt(re * re + im * im);
  double theta = acos(re / r);
  if (im < 0) theta = -theta;
  gsl_sf_result real, imaginary;
  gsl_sf_complex_dilog_e(r, theta, &real, &imaginary);
  return imaginary.val;
}

double log_re(const double re, const double im){//real part of the logrithm function
  gsl_sf_result real, imaginary;
  gsl_sf_complex_log_e(re, im, &real, &imaginary);
  return real.val;
}

double log_im(const double re, const double im){//imaginary part of the logrithm function
  gsl_sf_result real, imaginary;
  gsl_sf_complex_log_e(re, im, &real, &imaginary);
  return imaginary.val;
}

double f1(const double r){
  double t = sqrt(4.0 * r - 1.0);
  double result = (log(r) * log_im(1.0 - 1.0 / (2.0 * r), t / (2.0 * r)) + dilog_im(1.0 / (2.0 * r), t / (2.0 * r)) - dilog_im(1.0 / (2.0 * r), -t / (2.0 * r))) / t;
  return result;
}

double f2(const double r, const double r1, const double r2){
  double R = sqrt(r1 * r2);
  double rho = r1 / r2;
  if (rho == 1.0)
    return 0.5 * log(R) / R + (1.0 - r + 0.5 * r * log(r)) / (1.0 - r) / R;
  double a = sqrt(rho) * log(rho) / (2.0 * (rho - 1.0));
  double b = sqrt(rho) / (2.0 * (rho - 1.0)) * ( r * log(r) * log(rho) / (1.0 - r) - dilog_re(1.0 - rho, 0.0) + dilog_re(1.0 - 1.0 / rho, 0.0));
  return a * log(R) / R + b / R;
}

int CharginoSVD(const double M2, const double mu, const double gu, const double gd){
  X(0,0) = M2; X(0,1) = gu * VEV / sqrt(2.0);
  X(1,0) = gd * VEV / sqrt(2.0); X(1,1) = mu;
  TDecompSVD svd(X);
  svd.Decompose();
  U = svd.GetU();
  V = svd.GetV();
  MX = svd.GetSig();
  U.Invert();
  V.Invert();
  gDR(0,0) = gu * V(0,1) * U(0,0) + gd * V(0,0) * U(0,1);
  gDR(0,1) = gu * V(0,1) * U(1,0) + gd * V(0,0) * U(1,1);
  gDR(1,0) = gu * V(1,1) * U(0,0) + gd * V(1,0) * U(0,1);
  gDR(1,1) = gu * V(1,1) * U(1,0) + gd * V(1,0) * U(1,1);
  gDL(0,0) = gDR(0,0);
  gDL(0,1) = gDR(1,0);
  gDL(1,0) = gDR(0,1);
  gDL(1,1) = gDR(1,1);
  double s2 = 0.5 - pow(sin(thetaW), 2);
  double c2 = pow(cos(thetaW), 2);
  GR(0,0) = c2 * U(0,0) * U(0,0) + s2 * U(0,1) * U(0,1);
  GR(0,1) = c2 * U(0,0) * U(1,0) + s2 * U(0,1) * U(1,1);
  GR(1,0) = c2 * U(1,0) * U(0,0) + s2 * U(1,1) * U(0,1);
  GR(1,1) = c2 * U(1,0) * U(1,0) + s2 * U(1,1) * U(1,1);
  GL(0,0) = c2 * V(0,0) * V(0,0) + s2 * V(0,1) * V(0,1);
  GL(0,1) = c2 * V(0,0) * V(1,0) + s2 * V(0,1) * V(1,1);
  GL(1,0) = c2 * V(1,0) * V(0,0) + s2 * V(1,1) * V(0,1);
  GL(1,1) = c2 * V(1,0) * V(1,0) + s2 * V(1,1) * V(1,1);
  return 0;
}

double d1(const double M2, const double mu, const double gu, const double gd){
  CharginoSVD(M2, mu, gu, gd);
  double r0 = pow(MX(0) / MH, 2);
  double r1 = pow(MX(1) / MH, 2);
  double factor = pow(alpha, 2) / (4.0 * sqrt(2.0) * M_PI * M_PI * pow(sin(thetaW), 2)) / (MW * pow(MH, 2));
  double g = 2.0 * MW / VEV;
  double eta0 = pow(pdf->alphasQ(MX(0)) / alphaS0, 8.0 / 46.0);
  double eta1 = pow(pdf->alphasQ(MX(1)) / alphaS0, 8.0 / 46.0);
  double result = factor * (gDR(0,0) * MX(0) * f1(r0) * pow(eta0, 4) + gDR(1,1) * MX(1) * f1(r1) * pow(eta1, 4)) / g;
  return result;
}

double d1F(const double M2, const double mu, const double gu, const double gd){
  CharginoSVD(M2, mu, gu, gd);
  double r0 = pow(MX(0) / MH, 2);
  double r1 = pow(MX(1) / MH, 2);
  double factor = alpha / (4.0 * 2.0 * M_PI * M_PI * M_PI) / (M2 * mu);
  double eta0 = pow(pdf->alphasQ(MX(0)) / alphaS0, 8.0 / 46.0);
  double eta1 = pow(pdf->alphasQ(MX(1)) / alphaS0, 8.0 / 46.0);
  double eta = sqrt(eta0 * eta1);
  double R = sqrt(r0 * r1);
  double rho = r0 / r1;
  double F1 = 0.0;
  if (rho == 1.0) F1 = -0.5 * log(R) - 1.0 + 0.5;
  else F1 = -0.5 * log(R) - 1.0 + (rho + 1.0) * log(rho) / (rho - 1.0) / 4.0;
  double result = factor * gu * gd * F1 * pow(eta, 4);
  return result;
}

double d2(const double M2, const double mu, const double gu, const double gd){
  CharginoSVD(M2, mu, gu, gd);
  double r = pow(MZ / MH, 2);
  double r0 = pow(MX(0) / MH, 2);
  double r1 = pow(MX(1) / MH, 2);
  double g = 2.0 * MW / VEV;
  double eta0 = pow(pdf->alphasQ(MX(0)) / alphaS0, 8.0 / 46.0);
  double eta1 = pow(pdf->alphasQ(MX(1)) / alphaS0, 8.0 / 46.0);
  double factor = pow(alpha, 2) / (16.0 * sqrt(2.0) * M_PI * M_PI * pow(cos(thetaW) * sin(thetaW) * sin(thetaW), 2)) / (MW * pow(MH, 2));
  double term0 = (gDR(0,0) * GR(0,0)) * f2(r, r0, r0) + (gDR(0,1) * GR(1,0)) * f2(r, r0, r1);
  double term1 = (gDR(1,0) * GR(0,1)) * f2(r, r1, r0) + (gDR(1,1) * GR(1,1)) * f2(r, r1, r1);
  double result = factor * (term0 * MX(0) * pow(eta0, 4) + term1 * MX(1) * pow(eta1, 4)) / g;
  return sqrt(2.0) * result;
}

double d2F(const double M2, const double mu, const double gu, const double gd){
  CharginoSVD(M2, mu, gu, gd);
  double r = pow(MZ / MH, 2);
  double r0 = pow(MX(0) / MH, 2);
  double r1 = pow(MX(1) / MH, 2);
  double factor = alpha / (16.0 * M_PI * M_PI * M_PI * pow(cos(thetaW), 2)) / (M2 * mu);
  double eta0 = pow(pdf->alphasQ(MX(0)) / alphaS0, 8.0 / 46.0);
  double eta1 = pow(pdf->alphasQ(MX(1)) / alphaS0, 8.0 / 46.0);
  double eta = sqrt(eta0 * eta1);
  double R = sqrt(r0 * r1);
  double rho = r0 / r1;
  double A2 = 0.0;
  double B2 = 0.0;
  if (rho == 1.0){
    A2 = -3.0 / 8.0 + 0.5 * pow(sin(thetaW), 2);
    B2 = - (4.0 * pow(sin(thetaW), 2) - 3.0) * (1.0 - r + r * log(r)) / (8.0 * (r - 1.0));
  }
  else {
    A2 = ((rho - 1.0) * (2.0 - rho) - rho * log(rho)) / (4.0 * pow(rho - 1.0, 2)) + 0.5 * pow(sin(thetaW), 2);
    B2 = ( (2.0 - 2.0 * r + r * log(r)) * (rho - 1.0) * (rho - 2.0 - 2.0 * pow(sin(thetaW), 2) * (rho - 1.0)) + (r - 1.0) * (rho - 1.0) * (0.5 * rho + 1.0 - pow(sin(thetaW), 2) * (rho + 1.0)) * log(rho) + r * rho * log(r) * log(rho) + (r - 1.0) * rho * (dilog_re(1.0 - rho, 0.0) - dilog_re(1.0 - 1.0 / rho, 0.0))) / (4.0 * (r - 1.0) * pow(rho - 1.0, 2));
  }
  double F2 = A2 * log(R) + B2;
  double result = factor * gu * gd * F2 * pow(eta, 4);
  return result;
}

double d1u(const double M2, const double mu, const double gu, const double gd){
  double factor = d1F(M2, mu, gu, gd);
  return factor * (2.0 / 3.0) * Mu * 0.197e-13;//in unit e.cm
}

double d1d(const double M2, const double mu, const double gu, const double gd){
  double factor = d1F(M2, mu, gu, gd);
  return factor * (-1.0 / 3.0) * Md * 0.197e-13;//in unit e.cm
}

double d2u(const double M2, const double mu, const double gu, const double gd){
  double factor = d2F(M2, mu, gu, gd);
  return factor * (0.5 - 2.0 * pow(sin(thetaW), 2) * 2.0 / 3.0) * Mu * 0.197e-13;//in unit e.cm
}

double d2d(const double M2, const double mu, const double gu, const double gd){
  double factor = d2F(M2, mu, gu, gd);
  return factor * (-0.5 + 2.0 * pow(sin(thetaW), 2) / 3.0) * Md * 0.197e-13;//in unit e.cm
}

double dp(const double du, const double dd, const double * gt){
  double gtu = gt[0];
  double gtd = gt[1];
  double result = du * gtu + dd * gtd;
  return result;
}

double dperror(const double du, const double dd, const double * et){
  double etu = et[0];
  double etd = et[1];
  double etud = et[2];
  double err2 = pow(du * etu, 2) + pow(dd * etd, 2) + 2.0 * du * dd * etud;
  return sqrt(err2);
}

double dn(const double du, const double dd, const double * gt){
  double gtu = gt[0];
  double gtd = gt[1];
  double result = du * gtd + dd * gtu;
  return result;
}

double dnerror(const double du, const double dd, const double * et){
  double etu = et[0];
  double etd = et[1];
  double etud = et[2];
  double err2 = pow(du * etd, 2) + pow(dd * etu, 2) + 2.0 * du * dd * etud;
  return sqrt(err2);
}

double nEDMlimit = 3.0e-28 / 2.0;
double Solve_mu_nEDM(const double M2, const double * par){
  double gu = par[0];
  double gd = par[1];
  double gt[2] = {par[2], par[3]};
  double et[3] = {par[4], par[5], par[6]};
  double Min = log10(190.0);
  double Max = log10(2.0e4);
  double mass = Min;
  double du = d1u(M2, pow(10.0, mass), gu, gd) + d2u(M2, pow(10.0, mass), gu, gd);//
  double dd = d1d(M2, pow(10.0, mass), gu, gd) + d2d(M2, pow(10.0, mass), gu, gd);//
  double nedm = abs(dn(du, dd, gt));
  if (nedm < nEDMlimit) return pow(10.0, mass);
  double nedmerror = dnerror(du, dd, et);
  if (nedm - nedmerror < nEDMlimit) return pow(10.0, mass);
  while (true) {
    mass = mass + 0.001;
    du = d1u(M2, pow(10.0, mass), gu, gd);//
    dd = d1d(M2, pow(10.0, mass), gu, gd);//
    nedm = abs(dn(du, dd, gt));
    nedmerror = dnerror(du, dd, et);
    //cout << "...  " << nedm << " " << nedmerror << endl;
    if (nedm - nedmerror < nEDMlimit || mass > Max) break;
  }
  return pow(10.0, mass);
}

int Test(){
  TMatrixD X(2,2);
  X(0,0) = 1.0; X(0,1) = 2.0;
  X(1,0) = 3.0; X(1,1) = 4.0;
  TDecompSVD Xsvd(X);
  Xsvd.Decompose();
  TMatrixD U = Xsvd.GetU();
  TMatrixD V = Xsvd.GetV();
  TVectorD M = Xsvd.GetSig();
  U.Print();
  M.Print();
  V.Print();
  return 0;
}





#endif
