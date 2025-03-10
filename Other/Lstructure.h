#ifndef _LSTRUCTURE_H_
#define _LSTRUCTURE_H_

#include <iostream>
#include <cstring>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

class Lstructure{
 public:
  int f1_ct10(double * var, double * f1);/* unpolarized PDF */
  int f1_mmht(double * var, double * f1);
  int f1p(double * var, double * f1, const char set[]);
  int f1n(double * var, double * f1, const char set[]);
  int f1N(double * AZ, double * var, double * f1, const char set[]);
  int g1_grsv(double * var, double * g1);/* helicity */
  int g1_dssv(double * var, double * g1);
  int g1p(double * var, double * g1, const char set[]);
  int g1n(double * var, double * g1, const char set[]);
  int g1N(double * AZ, double * var, double * g1, const char set[]);
  int h1_std(double * var, double * h1, double * para);/* transversity */
  int h1p(double * var, double * h1, double * para);
  int h1n(double * var, double * h1, double * para);
  int h1N(double * AZ, double * var, double * h1, double * para);
  int f1t_std(double * var, double * f1t, double * para);/* Sivers */
  int f1tM_std(double * var, double * f1tM, double * para);
  int f1tp(double * var, double * f1t, double * para);
  int f1tn(double * var, double * f1t, double * para);
  int f1tN(double * AZ, double * var, double * f1t, double * para);
  int h1t_std(double * var, double * h1t, double * para);/* pretzelosity */
  int h1tM_std(double * var, double * h1tM, double * para);
  int h1tp(double * var, double * h1t, double * para);
  int h1tn(double * var, double * h1t, double * para);
  int h1tN(double * AZ, double * var, double * h1t, double * para);
  int D1_dss(double * var, double * D1);
  int D1pip(double * var, double * D1, const char set[]);
  int D1pim(double * var, double * D1, const char set[]);
  int H1_std(double * var, double * H1, double * para);
  int H1_poly(double * var, double * H1, double * para);
  int H1pip(double * var, double * H1, double * para, const char set[]);
  int H1pim(double * var, double * H1, double * para, const char set[]);
  int FUUTp(double * var, double * FF, double * fpara, double * Dpara);
  int FUUTn(double * var, double * FF, double * fpara, double * Dpara);
  int FUUTN(double * AZ, double * var, double * FF, double * fpara, double * Dpara);
  int FUTsinHmSp(double * var, double * FF, double * fpara, double * Dpara);
  int FUTsinHmSn(double * var, double * FF, double * fpara, double * Dpara);
  int FUTsinHmSN(double * AZ, double * var, double * FF, double * fpara, double * Dpara);
  int FUTsinHpSp(double * var, double * FF, double * hpara, double * Hpara);
  int FUTsinHpSn(double * var, double * FF, double * hpara, double * Hpara);
  int FUTsinHpSN(double * AZ, double * var, double * FF, double * hpara, double * Hpara);
  int FUTsin3HmSp(double * var, double * FF, double * hpara, double * Hpara);
  int FUTsin3HmSn(double * var, double * FF, double * hpara, double * Hpara);
  int FUTsin3HmSN(double * AZ, double * var, double * FF, double * hpara, double * Hpara);
  int AsinHmSp(double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int AsinHmSn(double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int AsinHmSN(double * AZ, double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int AsinHpSp(double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int AsinHpSn(double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int AsinHpSN(double * AZ, double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int Asin3HmSp(double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int Asin3HmSn(double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  int Asin3HmSN(double * AZ, double * kin, double * Asym, double * Tf, double * TD, double * Uf, double * UD);
  double xsprefactor(double * kin);
  int CalcVariables(double * lab, double * phys);
  double Jacobian(double * lab);//d trento / d lab
  int sigmaUUT(double * AZ, double * lab, double * xs, double * fpara, double * Dpara);// lab frame
};

Lstructure * lst;

const LHAPDF::PDF * _pdfs = LHAPDF::mkPDF("CT10", 0);
int Lstructure::f1_ct10(double * var, double * f1){
  //var: x, Q2
  double x = var[0];
  double Q2 = var[1];
  //f1: u, d, s, ubar, dbar, sbar
  f1[0] = _pdfs->xfxQ2(2, x, Q2) / x;//u
  f1[1] = _pdfs->xfxQ2(1, x, Q2) / x;//d
  f1[2] = _pdfs->xfxQ2(3, x, Q2) / x;//s
  f1[3] = _pdfs->xfxQ2(-2, x, Q2) / x;//ubar
  f1[4] = _pdfs->xfxQ2(-1, x, Q2) / x;//dbar
  f1[5] = _pdfs->xfxQ2(-3, x, Q2) / x;//sbar
  return 0;
}

const LHAPDF::PDF* _pdfMMHT = LHAPDF::mkPDF("MMHT2014lo68cl", 0);
int Lstructure::f1_mmht(double * var, double * f1){
  //var: x, Q2
  double x = var[0];
  double Q2 = var[1];
  //f1: u, d, s, ubar, dbar, sbar
  f1[0] = _pdfMMHT->xfxQ2(2, x, Q2) / x;//u
  f1[1] = _pdfMMHT->xfxQ2(1, x, Q2) / x;//d
  f1[2] = _pdfMMHT->xfxQ2(3, x, Q2) / x;//s
  f1[3] = _pdfMMHT->xfxQ2(-2, x, Q2) / x;//ubar
  f1[4] = _pdfMMHT->xfxQ2(-1, x, Q2) / x;//dbar
  f1[5] = _pdfMMHT->xfxQ2(-3, x, Q2) / x;//sbar
  return 0;
}

int Lstructure::f1p(double * var, double * f1, const char set[] = "ct10"){
  if (strcmp(set,"ct10") == 0)
    lst->f1_ct10(var, f1);
  else if (strcmp(set,"mmht") == 0)
    lst->f1_mmht(var, f1);
  else {
    std::cout << "Lstructure::f1p: No such PDF set!" << std::endl;
    return 1;
  }
  return 0;
}

int Lstructure::f1n(double * var, double * f1, const char set[] = "ct10"){
  double pf1[6];
  if (strcmp(set,"ct10") == 0)
    lst->f1_ct10(var, pf1);
  else if (strcmp(set,"mmht") == 0)
    lst->f1_mmht(var, pf1);
  else {
    std::cout << "Lstructure::f1n: No such PDF set!" << std::endl;
    return 1;
  }
  f1[0] = pf1[1];
  f1[1] = pf1[0];
  f1[2] = pf1[2];
  f1[3] = pf1[4];
  f1[4] = pf1[3];
  f1[5] = pf1[5];
  return 0;
}

int Lstructure::f1N(double * AZ, double * var, double * f1, const char set[] = "ct10"){
  //AZ: Np, Nn
  double Np = AZ[0];
  double Nn = AZ[1];
  double pf1[6], nf1[6];
  lst->f1p(var, pf1, set);
  lst->f1n(var, nf1, set);
  for (int i = 0; i < 6; i++){
    f1[i] = Np * pf1[i] + Nn * nf1[i];
  }
  return 0;
}


extern "C"{
  extern struct { double IINI;} gintini_;
}
extern "C"{
  void parpol_(int * ISET, double * X, double * Q2, double * U, double * D, double * UB, double * DB, double * ST, double * GL, double * G1P, double * G1N);
}
int Lstructure::g1_grsv(double * var, double * g1){
  //var: x, Q2
  double x = var[0];
  double Q2 = var[1];
  int iset = 3;
  double value[8];
  parpol_(&iset, &x, &Q2, &value[0], &value[1], &value[2], &value[3], &value[4], &value[5], &value[6], &value[7]);
  g1[0] = value[0] / x;//u
  g1[1] = value[1] / x;//d
  g1[2] = value[4] / 2.0 / x;//s
  g1[3] = value[2] / x;//ubar
  g1[4] = value[3] / x;//dbar
  g1[5] = value[4] / 2.0 / x;//sbar
  return 0;
}

extern "C"{
  void dssvini_(int * iset);
  void dssvfit_(double * X, double * Q2, double * xDUV, double * xDDV, double * xDUBAR, double * xDDBAR, double * xDSTR, double * xDGLU);
}
bool dssvinit = true;
int Lstructure::g1_dssv(double * var, double * g1){
  //var: x, Q2
  double x = var[0];
  double Q2 = var[1];
  int iset;
  if (dssvinit){
    iset = 0;
    dssvini_(&iset);
    dssvinit = false;
  }
  double value[6];
  dssvfit_(&x, &Q2, &value[0], &value[1], &value[2], &value[3], &value[4], &value[5]);
  g1[0] = (value[0] + value[2]) / x;//u
  g1[1] = (value[1] + value[3]) / x;//d
  g1[2] = value[4] / x;//s
  g1[3] = value[2] / x;//ubar
  g1[4] = value[3] / x;//dbar
  g1[5] = value[4] / x;//sbar
  return 0;
}

int Lstructure::g1p(double * var, double * g1, const char set[] = "grsv"){
  if (strcmp(set,"grsv") == 0)
    lst->g1_grsv(var, g1);
  else if (strcmp(set,"dssv") == 0)
    lst->g1_dssv(var, g1);
  else {
    std::cout << "Lstructure::g1p: No such g1 set!" << std::endl;
    return 1;
  }
  return 0;
}

int Lstructure::g1n(double * var, double * g1, const char set[] = "grsv"){
  double pg1[6];
  if (strcmp(set,"grsv") == 0)
    lst->g1_grsv(var, pg1);
  else if (strcmp(set,"dssv") == 0)
    lst->g1_dssv(var, pg1);
  else {
    std::cout << "Lstructure::g1n: No such g1 set!" << std::endl;
    return 1;
  }
  g1[0] = pg1[1];
  g1[1] = pg1[0];
  g1[2] = pg1[2];
  g1[3] = pg1[4];
  g1[4] = pg1[3];
  g1[5] = pg1[5];
  return 0;
}

int Lstructure::g1N(double * AZ, double * var, double * g1, const char set[] = "grsv"){
  //AZ: Np, Nn, Polp, Poln
  double Np = AZ[0];
  double Nn = AZ[1];
  double Polp = AZ[2];
  double Poln = AZ[3];
  double pg1[6];
  double ng1[6];
  lst->g1p(var, pg1, set);
  lst->g1n(var, ng1, set);
  for (int i = 0; i < 6; i++){
    g1[i] = Np * Polp * pg1[i] + Nn * Poln * ng1[i];
  }
  return 0;
}

int Lstructure::h1_std(double * var, double * h1, double * para = 0){
  //"standard" parametrization, Anselmino et al. PR D87 (2013) 094019
  double Nu, Nd, au, ad, bu, bd, kt2;
  if (para == 0){
    Nu = 0.46;
    Nd = -1.0;
    au = 1.11;
    ad = 1.11;
    bu = 3.64;
    bd = 3.64;
    kt2 = 0.25;
  }
  else {
    Nu = para[0];
    Nd = para[1];
    au = para[2];
    ad = para[3];
    bu = para[4];
    bd = para[5];
    kt2 = para[6];
  }
  double x = var[0];
  double f1[6], g1[6];
  lst->f1p(var, f1, "ct10");
  lst->g1p(var, g1, "grsv");
  h1[0] = 0.5 * Nu * pow(x,au) * pow(1.0-x,bu) * pow(au+bu,au+bu)/pow(au,au)/pow(bu,bu) * (f1[0] + g1[0]);//u
  h1[1] = 0.5 * Nd * pow(x,ad) * pow(1.0-x,bd) * pow(ad+bd,ad+bd)/pow(ad,ad)/pow(bd,bd) * (f1[1] + g1[1]);//d
  h1[2] = 0;
  h1[3] = 0;
  h1[4] = 0;
  h1[5] = 0;
  return 0;
}

int Lstructure::h1p(double * var, double * h1, double * para = 0){
  lst->h1_std(var, h1, para);
  return 0;
}
int Lstructure::h1n(double * var, double * h1, double * para = 0){
  double ph1[6];
  lst->h1_std(var, ph1, para);
  h1[0] = ph1[1];
  h1[1] = ph1[0];
  h1[2] = ph1[2];
  h1[3] = ph1[4];
  h1[4] = ph1[3];
  h1[5] = ph1[5];
  return 0;
}

int Lstructure::h1N(double * AZ, double * var, double * h1, double * para = 0){
  //AZ: Np, Nn, Polp, Poln
  double Np = AZ[0];
  double Nn = AZ[1];
  double Polp = AZ[2];
  double Poln = AZ[3];
  double ph1[6], nh1[6];
  lst->h1p(var, ph1, para);
  lst->h1n(var, nh1, para);
  for (int i = 0; i < 6; i++){
    h1[i] = Np * Polp * ph1[i] + Nn * Poln * nh1[i];
  }
  return 0;
}

int Lstructure::f1t_std(double * var, double * f1t, double * para = 0){
  //var: x, Q2
  double x = var[0];
  double Nu, Nd, Ns, Nubar, Ndbar, Nsbar;
  double au, ad, asea, b, Ms2;
  if (para == 0){
    Nu = 0.35;
    Nd = -0.90;
    Ns = -0.24;
    Nubar = 0.04;
    Ndbar = -0.40;
    Nsbar = 1.0;
    au = 0.73;
    ad = 1.08;
    asea = 0.79;
    b = 3.46;
    Ms2 = 0.34;
  }
  else {
    Nu = para[0];
    Nd = para[1];
    Ns = para[2];
    Nubar = para[3];
    Ndbar = para[4];
    Nsbar = para[5];
    au = para[6];
    ad = para[7];
    asea = para[8];
    b = para[9];
    Ms2 = para[10];
  }
  const double Mproton = 0.93827;
  double f1[6];
  lst->f1p(var, f1, "ct10");
  f1t[0] = - f1[0] * Nu * pow(x,au) * pow(1.0-x,b) * pow(au+b,au+b) / pow(au,au) / pow(b,b) * sqrt(2.0*exp(1)) * Mproton * sqrt(Ms2) / (Ms2 + 0.25);//u
  f1t[1] = - f1[1] * Nd * pow(x,ad) * pow(1.0-x,b) * pow(ad+b,ad+b) / pow(ad,ad) / pow(b,b) * sqrt(2.0*exp(1)) * Mproton * sqrt(Ms2) / (Ms2 + 0.25);//d
  f1t[2] = - f1[2] * Ns * pow(x,asea) * pow(1.0-x,b) * pow(asea+b,asea+b) / pow(asea,asea) / pow(b,b) * sqrt(2.0*exp(1)) * Mproton * sqrt(Ms2) / (Ms2 + 0.25);//s
  f1t[3] = - f1[3] * Nubar * pow(x,asea) * pow(1.0-x,b) * pow(asea+b,asea+b) / pow(asea,asea) / pow(b,b) * sqrt(2.0*exp(1)) * Mproton * sqrt(Ms2) / (Ms2 + 0.25);//ubar
  f1t[4] = - f1[4] * Ndbar * pow(x,asea) * pow(1.0-x,b) * pow(asea+b,asea+b) / pow(asea,asea) / pow(b,b) * sqrt(2.0*exp(1)) * Mproton * sqrt(Ms2) / (Ms2 + 0.25);//dbar
  f1t[5] = - f1[5] * Nsbar * pow(x,asea) * pow(1.0-x,b) * pow(asea+b,asea+b) / pow(asea,asea) / pow(b,b) * sqrt(2.0*exp(1)) * Mproton * sqrt(Ms2) / (Ms2 + 0.25);//sbar
  return 0;
}

int Lstructure::f1tM_std(double * var, double * f1tM, double * para = 0){
  double Ms2;
  if (para == 0) Ms2 = 0.34;
  else Ms2 = para[10];
  double kt2 = 0.25 * Ms2 / (Ms2 + 0.25);
  const double Mproton = 0.93827;
  double f1t[6];
  lst->f1t_std(var, f1t, para);
  for (int i = 0; i < 6; i++){
    f1tM[i] = kt2 / (2.0 * pow(Mproton,2)) * f1t[i];
  }
  return 0;
}

int Lstructure::f1tp(double * var, double * f1t, double * para = 0){
  lst->f1t_std(var, f1t, para);
  return 0;
}

int Lstructure::f1tn(double * var, double * f1t, double * para = 0){
  double pf1t[6];
  lst->f1t_std(var, pf1t, para);
  f1t[0] = pf1t[1];
  f1t[1] = pf1t[0];
  f1t[2] = pf1t[2];
  f1t[3] = pf1t[4];
  f1t[4] = pf1t[3];
  f1t[5] = pf1t[5];
  return 0;
}

int Lstructure::f1tN(double * AZ, double * var, double * f1t, double * para){
  //AZ: Np, Nn, Polp, Poln
  double Np = AZ[0];
  double Nn = AZ[1];
  double Polp = AZ[2];
  double Poln = AZ[3];
  double pf1t[6], nf1t[6];
  lst->f1tp(var, pf1t, para);
  lst->f1tn(var, nf1t, para);
  for (int i = 0; i < 6; i++){
    f1t[i] = Np * Polp * pf1t[i] + Nn * Poln * nf1t[i];
  }
  return 0;
}

int Lstructure::h1t_std(double * var, double * h1t, double * para = 0){
  //var: x, Q2
  double x = var[0];
  double Nu, Nd, au, ad, bu, bd, Mt2;
  if (para == 0){
    Nu = 1.0;
    Nd = -1.0;
    au = 2.5;
    ad = 2.5;
    bu = 2.0;
    bd = 2.0;
    Mt2 = 0.18;
  }
  else {
    Nu = para[0];
    Nd = para[1];
    au = para[2];
    ad = para[3];
    bu = para[4];
    bd = para[5];
    Mt2 = para[6];
  }
  double f1[6], g1[6];
  const double Mproton = 0.93827;
  lst->f1p(var, f1, "ct10");
  lst->g1p(var, g1, "dssv");
  h1t[0] = (f1[0] - g1[0]) * Nu * pow(x,au) * pow(1.0-x,bu) * pow(au+bu,au+bu) / pow(au,au) / pow(bu,bu) * exp(1) * pow(Mproton,2) / (Mt2 + 0.25);//u
  h1t[1] = (f1[1] - g1[1]) * Nd * pow(x,ad) * pow(1.0-x,bd) * pow(ad+bd,ad+bd) / pow(ad,ad) / pow(bd,bd) * exp(1) * pow(Mproton,2) / (Mt2 + 0.25);//d
  h1t[2] = 0;
  h1t[3] = 0;
  h1t[4] = 0;
  h1t[5] = 0;
  return 0;
}

int Lstructure::h1tM_std(double * var, double * h1tM, double * para = 0){
  double Mt2;
  if (para == 0) Mt2 = 0.18;
  else Mt2 = para[6];
  double kt2 = 0.25 * Mt2 / (Mt2 + 0.25);
  double h1t[6];
  const double Mproton = 0.93827;
  lst->h1t_std(var, h1t, para);
  for (int i = 0; i < 6; i++){
    h1tM[i] = kt2 / (2.0 * pow(Mproton,2)) * h1t[i];
  }
  return 0;
}

int Lstructure::h1tp(double * var, double * h1t, double * para = 0){
  lst->h1t_std(var, h1t, para);
  return 0;
}

int Lstructure::h1tn(double * var, double * h1t, double * para = 0){
  double ph1t[6];
  lst->h1t_std(var, ph1t, para);
  h1t[0] = ph1t[1];
  h1t[1] = ph1t[0];
  h1t[2] = ph1t[2];
  h1t[3] = ph1t[4];
  h1t[4] = ph1t[3];
  h1t[5] = ph1t[5];
  return 0;
}

int Lstructure::h1tN(double * AZ, double * var, double * h1t, double * para = 0){
  //AZ: Np, Nn, Polp, Poln
  double Np = AZ[0];
  double Nn = AZ[1];
  double Polp = AZ[2];
  double Poln = AZ[3];
  double ph1t[6], nh1t[6];
  lst->h1tp(var, ph1t, para);
  lst->h1tn(var, nh1t, para);
  for (int i = 0; i < 6; i++){
    h1t[i] = Np * Polp * ph1t[i] + Nn * Poln * nh1t[i];
  }
  return 0;
}

extern "C"{
  extern struct { int FINI;} fraginid_;
}
extern "C"{
  void fdss_(int * IH, int * IC, int * IO, double * Z, double * Q2, double * zU, double * zUBAR, double * zD, double * zDBAR, double * zS, double * zSBAR, double * zC, double * zB, double * zG);
}
int Lstructure::D1_dss(double * var, double * D1){
  //var: z, Q2
  double z = var[0];
  double Q2 = var[1];
  int ih = 1;
  int ic = 1;
  int io = 0;
  double u, ubar, d, dbar, s, sbar, c, b, g;
  fdss_ (&ih, &ic, &io, &z, &Q2, &u, &ubar, &d, &dbar, &s, &sbar, &c, &b, &g);
  D1[0] = u / z;//u
  D1[1] = d / z;//d
  D1[2] = s / z;//s
  D1[3] = ubar / z;//ubar
  D1[4] = dbar / z;//dbar
  D1[5] = sbar / z;//sbar
  return 0;
}

int Lstructure::D1pip(double * var, double * D1, const char set[] = "dss"){
  if (strcmp(set,"dss") == 0)
    lst->D1_dss(var, D1);
  else {
    std::cout << "Lstructure::D1pip: No such D1 set!" << std::endl;
    return 1;
  }
  return 0;
}

int Lstructure::D1pim(double * var, double * D1, const char set[] = "dss"){
  double D1p[6];
  if (strcmp(set,"dss") == 0)
    lst->D1_dss(var, D1p);
  else {
    std::cout << "Lstructure::D1pim: No such D1 set!" << std::endl;
    return 1;
  }
  D1[0] = D1p[3];
  D1[1] = D1p[4];
  D1[2] = D1p[5];
  D1[3] = D1p[0];
  D1[4] = D1p[1];
  D1[5] = D1p[2];
  return 0;
}

int Lstructure::H1_std(double * var, double * H1, double * para = 0){
  //var: z, Q2
  double z = var[0];
  double Nu, Nd, cu, cd, du, dd, Mh2;
  if (para == 0){
    Nu = 0.49;
    Nd = -1.0;
    cu = 1.06;
    cd = 1.06;
    du = 0.07;
    dd = 0.07;
    Mh2 = 1.5;
  }
  else {
    Nu = para[0];
    Nd = para[1];
    cu = para[2];
    cd = para[3];
    du = para[4];
    dd = para[5];
    Mh2 = para[6];
  }
  const double Mpion = 0.13957;
  double D1[6];
  lst->D1_dss(var, D1);
  H1[0] = D1[0] * Nu * z * pow(z,cu) * pow(1.0-z,du) * pow(cu+du,cu+du) / pow(cu,cu) / pow(du,du) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//u
  H1[1] = D1[1] * Nd * z * pow(z,cd) * pow(1.0-z,dd) * pow(cd+dd,cd+dd) / pow(cd,cd) / pow(dd,dd) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//d
  H1[2] = D1[2] * Nd * z * pow(z,cd) * pow(1.0-z,dd) * pow(cd+dd,cd+dd) / pow(cd,cd) / pow(dd,dd) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//s
  H1[3] = D1[3] * Nd * z * pow(z,cd) * pow(1.0-z,dd) * pow(cd+dd,cd+dd) / pow(cd,cd) / pow(dd,dd) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//ubar
  H1[4] = D1[4] * Nu * z * pow(z,cu) * pow(1.0-z,du) * pow(cu+du,cu+du) / pow(cu,cu) / pow(du,du) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//dbar
  H1[5] = D1[5] * Nd * z * pow(z,cd) * pow(1.0-z,dd) * pow(cd+dd,cd+dd) / pow(cd,cd) / pow(dd,dd) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//sbar
  return 0;
}

int Lstructure::H1_poly(double * var, double * H1, double * para = 0){
  //var: z, Q2
  double z = var[0];
  double Nu, Nd, cu, cd, du, dd, Mh2;
  if (para == 0){
    Nu = 1.0;
    Nd = -1.0;
    cu = -2.36;
    cd = -2.36;
    du = 2.12;
    dd = 2.12;
    Mh2 = 0.67;
  }
  else {
    Nu = para[0];
    Nd = para[1];
    cu = para[2];
    cd = para[3];
    du = para[4];
    dd = para[5];
    Mh2 = para[6];
  }
  const double Mpion = 0.13957;
  double D1[6];
  lst->D1_dss(var, D1);
  H1[0] = D1[0] * Nu * z * z * ((1.0 - cu - du) + cu * z + du * z * z) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//u
  H1[1] = D1[1] * Nd * z * z * ((1.0 - cd - dd) + cd * z + dd * z * z) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//d
  H1[2] = D1[2] * Nd * z * z * ((1.0 - cd - dd) + cd * z + dd * z * z) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//s
  H1[3] = D1[3] * Nd * z * z * ((1.0 - cd - dd) + cd * z + dd * z * z) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//ubar
  H1[4] = D1[4] * Nu * z * z * ((1.0 - cu - du) + cu * z + du * z * z) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//dbar
  H1[5] = D1[5] * Nd * z * z * ((1.0 - cd - dd) + cd * z + dd * z * z) * sqrt(2.0*exp(1)) * Mpion * sqrt(Mh2) / (Mh2 + 0.20);//sbar
  return 0;
}

int Lstructure::H1pip(double * var, double * H1, double * para = 0, const char set[] = "std"){
  if (strcmp(set, "std") == 0)
    lst->H1_std(var, H1, para);
  else if (strcmp(set, "poly") == 0)
    lst->H1_poly(var, H1, para);
  else {
    std::cout << "Lstructure::H1pip: No such H1 set!" << std::endl;
    return 1;
  }
  return 0;
}

int Lstructure::H1pim(double * var, double * H1, double * para = 0, const char set[] = "std"){
  double pH1[6];
  if (strcmp(set, "std") == 0)
    lst->H1_std(var, pH1, para);
  else if (strcmp(set, "poly") == 0)
    lst->H1_poly(var, pH1, para);
  else {
    std::cout << "Lstructure::H1pim: No such H1 set!" << std::endl;
    return 1;
  }
  H1[0] = pH1[3];
  H1[1] = pH1[4];
  H1[2] = pH1[5];
  H1[3] = pH1[0];
  H1[4] = pH1[1];
  H1[5] = pH1[2];
  return 0;
}

int Lstructure::FUUTp(double * var, double * FF, double * fpara = 0, double * Dpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double kt2, pt2;
  if (fpara == 0) kt2 = 0.25;
  else kt2 = fpara[0];
  if (Dpara == 0) pt2 = 0.20;
  else pt2 = Dpara[0];
  double Pt2 = pt2 + z*z*kt2;
  double f1[6], D1[6];
  double f1var[2] = {x, Q2};
  double D1var[2] = {z, Q2};
  lst->f1p(f1var, f1);
  lst->D1pip(D1var, D1);
  FF[0] = x * (pow(2.0/3.0, 2) * (f1[0]*D1[0] + f1[3]*D1[3]) + pow(1.0/3.0, 2) * (f1[1]*D1[1] + f1[2]*D1[2] + f1[4]*D1[4] + f1[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->D1pim(D1var, D1);
  FF[1] = x * (pow(2.0/3.0, 2) * (f1[0]*D1[0] + f1[3]*D1[3]) + pow(1.0/3.0, 2) * (f1[1]*D1[1] + f1[2]*D1[2] + f1[4]*D1[4] + f1[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}

int Lstructure::FUUTn(double * var, double * FF, double * fpara = 0, double * Dpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double kt2, pt2;
  if (fpara == 0) kt2 = 0.25;
  else kt2 = fpara[0];
  if (Dpara == 0) pt2 = 0.20;
  else pt2 = Dpara[0];
  double Pt2 = pt2 + z*z*kt2;
  double f1[6], D1[6];
  double f1var[2] = {x, Q2};
  double D1var[2] = {z, Q2};
  lst->f1n(f1var, f1);
  lst->D1pip(D1var, D1);
  FF[0] = x * (pow(2.0/3.0, 2) * (f1[0]*D1[0] + f1[3]*D1[3]) + pow(1.0/3.0, 2) * (f1[1]*D1[1] + f1[2]*D1[2] + f1[4]*D1[4] + f1[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->D1pim(D1var, D1);
  FF[1] = x * (pow(2.0/3.0, 2) * (f1[0]*D1[0] + f1[3]*D1[3]) + pow(1.0/3.0, 2) * (f1[1]*D1[1] + f1[2]*D1[2] + f1[4]*D1[4] + f1[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}

int Lstructure::FUUTN(double * AZ, double * var, double * FF, double * fpara = 0, double * Dpara = 0){
  //AZ: Np, Nn
  double Np = AZ[0];
  double Nn = AZ[1];
  double FFp[2], FFn[2];
  lst->FUUTp(var, FFp, fpara, Dpara);
  lst->FUUTn(var, FFn, fpara, Dpara);
  FF[0] = Np * FFp[0] + Nn * FFn[0];
  FF[1] = Np * FFp[1] + Nn * FFn[1];
  return 0;
}

int Lstructure::FUTsinHmSp(double * var, double * FF, double * fpara = 0, double * Dpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double Ms2, pt2;
  if (fpara == 0) Ms2 = 0.34;
  else Ms2 = fpara[10];
  if (Dpara == 0) pt2 = 0.20;
  else pt2 = Dpara[0];
  double kt2 = 0.25 * Ms2 / (Ms2 + 0.25);
  double Pt2 = pt2 + z*z*kt2;
  double f1t[6], D1[6];
  double fvar[2] = {x, Q2};
  double Dvar[2] = {z, Q2};
  lst->f1tp(fvar, f1t, fpara);
  lst->D1pip(Dvar, D1);
  const double Mproton = 0.93827;
  double pre = - z * kt2 * Pt / (Mproton * Pt2);
  FF[0] = pre * x * (pow(2.0/3.0,2) * (f1t[0]*D1[0] + f1t[3]*D1[3]) + pow(1.0/3.0,2) * (f1t[1]*D1[1] + f1t[2]*D1[2] + f1t[4]*D1[4] + f1t[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->D1pim(Dvar, D1);
  FF[1] = pre * x * (pow(2.0/3.0,2) * (f1t[0]*D1[0] + f1t[3]*D1[3]) + pow(1.0/3.0,2) * (f1t[1]*D1[1] + f1t[2]*D1[2] + f1t[4]*D1[4] + f1t[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}

int Lstructure::FUTsinHmSn(double * var, double * FF, double * fpara = 0, double * Dpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double Ms2, pt2;
  if (fpara == 0) Ms2 = 0.34;
  else Ms2 = fpara[10];
  if (Dpara == 0) pt2 = 0.20;
  else pt2 = Dpara[0];
  double kt2 = 0.25 * Ms2 / (Ms2 + 0.25);
  double Pt2 = pt2 + z*z*kt2;
  double f1t[6], D1[6];
  double fvar[2] = {x, Q2};
  double Dvar[2] = {z, Q2};
  lst->f1tn(fvar, f1t, fpara);
  lst->D1pip(Dvar, D1);
  const double Mneutron = 0.93957;
  double pre = - z * kt2 * Pt / (Mneutron * Pt2);
  FF[0] = pre * x * (pow(2.0/3.0,2) * (f1t[0]*D1[0] + f1t[3]*D1[3]) + pow(1.0/3.0,2) * (f1t[1]*D1[1] + f1t[2]*D1[2] + f1t[4]*D1[4] + f1t[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->D1pim(Dvar, D1);
  FF[1] = pre * x * (pow(2.0/3.0,2) * (f1t[0]*D1[0] + f1t[3]*D1[3]) + pow(1.0/3.0,2) * (f1t[1]*D1[1] + f1t[2]*D1[2] + f1t[4]*D1[4] + f1t[5]*D1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}

int Lstructure::FUTsinHmSN(double * AZ, double * var, double * FF, double * fpara = 0, double * Dpara = 0){
  //AZ: Np, Nn, Polp, Poln
  double Np = AZ[0];
  double Nn = AZ[1];
  double Polp = AZ[2];
  double Poln = AZ[3];
  double FFp[2], FFn[2];
  lst->FUTsinHmSp(var, FFp, fpara, Dpara);
  lst->FUTsinHmSn(var, FFn, fpara, Dpara);
  FF[0] = Np * Polp * FFp[0] + Nn * Poln * FFn[0];
  FF[1] = Np * Polp * FFp[1] + Nn * Poln * FFn[1];
  return 0;
}

int Lstructure::FUTsinHpSp(double * var, double * FF, double * hpara = 0, double * Hpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double hvar[2] = {x, Q2};
  double Hvar[2] = {z, Q2};
  double kt2, Mh2;
  if (hpara == 0) kt2 = 0.25;
  else kt2 = hpara[6];
  if (Hpara == 0) Mh2 = 0.67;
  else Mh2 = Hpara[6];
  double pt2 = 0.20 * Mh2 / (Mh2 + 0.20);
  double Pt2 = pt2 + z*z*kt2;
  double h1[6], H1[6];
  const double Mpion = 0.13957;
  double pre = Pt * pt2 / (z * Mpion * Pt2);
  lst->h1p(hvar, h1, hpara);
  lst->H1pip(Hvar, H1, Hpara);
  FF[0] = pre * x * (pow(2.0/3.0,2) * (h1[0]*H1[0] + h1[3]*H1[3]) + pow(1.0/3.0,2) * (h1[1]*H1[1] + h1[2]*H1[2] + h1[4]*H1[4] + h1[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->H1pim(Hvar, H1, Hpara);
  FF[1] = pre * x * (pow(2.0/3.0,2) * (h1[0]*H1[0] + h1[3]*H1[3]) + pow(1.0/3.0,2) * (h1[1]*H1[1] + h1[2]*H1[2] + h1[4]*H1[4] + h1[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}
 
int Lstructure::FUTsinHpSn(double * var, double * FF, double * hpara = 0, double * Hpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double hvar[2] = {x, Q2};
  double Hvar[2] = {z, Q2};
  double kt2, Mh2;
  if (hpara == 0) kt2 = 0.25;
  else kt2 = hpara[6];
  if (Hpara == 0) Mh2 = 0.67;
  else Mh2 = Hpara[6];
  double pt2 = 0.20 * Mh2 / (Mh2 + 0.20);
  double Pt2 = pt2 + z*z*kt2;
  double h1[6], H1[6];
  const double Mpion = 0.13957;
  double pre = Pt * pt2 / (z * Mpion * Pt2);
  lst->h1n(hvar, h1, hpara);
  lst->H1pip(Hvar, H1, Hpara);
  FF[0] = pre * x * (pow(2.0/3.0,2) * (h1[0]*H1[0] + h1[3]*H1[3]) + pow(1.0/3.0,2) * (h1[1]*H1[1] + h1[2]*H1[2] + h1[4]*H1[4] + h1[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->H1pim(Hvar, H1, Hpara);
  FF[1] = pre * x * (pow(2.0/3.0,2) * (h1[0]*H1[0] + h1[3]*H1[3]) + pow(1.0/3.0,2) * (h1[1]*H1[1] + h1[2]*H1[2] + h1[4]*H1[4] + h1[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}

int Lstructure::FUTsinHpSN(double * AZ, double * var, double * FF, double * hpara = 0, double * Hpara = 0){
  //AZ: Np, Nn, Polp, Poln
  double Np = AZ[0];
  double Nn = AZ[1];
  double Polp = AZ[2];
  double Poln = AZ[3];
  double FFp[2], FFn[2];
  lst->FUTsinHpSp(var, FFp, hpara, Hpara);
  lst->FUTsinHpSn(var, FFn, hpara, Hpara);
  FF[0] = Np * Polp * FFp[0] + Nn * Poln * FFn[0];
  FF[1] = Np * Polp * FFp[1] + Nn * Poln * FFn[1];
  return 0;
}

int Lstructure::FUTsin3HmSp(double * var, double * FF, double * hpara = 0, double * Hpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double hvar[2] = {x, Q2};
  double Hvar[2] = {z, Q2};
  double Mt2, Mh2;
  if (hpara == 0) Mt2 = 0.18;
  else Mt2 = hpara[6];
  if (Hpara == 0) Mh2 = 0.67;
  else Mh2 = Hpara[6];
  double kt2 = 0.25 * Mt2 / (Mt2 + 0.25);
  double pt2 = 0.20 * Mh2 / (Mh2 + 0.20);
  double Pt2 = pt2 + z*z*kt2;
  double h1t[6], H1[6];
  const double Mproton = 0.93827;
  const double Mpion = 0.13957;
  double pre = 0.5 * z * pow(kt2,2) * pt2 * pow(Pt,3) / (pow(Mproton,2) * Mpion * pow(Pt2,3));
  lst->h1tp(hvar, h1t, hpara);
  lst->H1pip(Hvar, H1, Hpara);
  FF[0] = pre * x * (pow(2.0/3.0,2) * (h1t[0]*H1[0] + h1t[3]*H1[3]) + pow(1.0/3.0,2) * (h1t[1]*H1[1] + h1t[2]*H1[2] + h1t[4]*H1[4] + h1t[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->H1pim(Hvar, H1, Hpara);
  FF[1] = pre * x * (pow(2.0/3.0,2) * (h1t[0]*H1[0] + h1t[3]*H1[3]) + pow(1.0/3.0,2) * (h1t[1]*H1[1] + h1t[2]*H1[2] + h1t[4]*H1[4] + h1t[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}

int Lstructure::FUTsin3HmSn(double * var, double * FF, double * hpara = 0, double * Hpara = 0){
  //var: x, Q2, z, Pt
  double x = var[0];
  double Q2 = var[1];
  double z = var[2];
  double Pt = var[3];
  double hvar[2] = {x, Q2};
  double Hvar[2] = {z, Q2};
  double Mt2, Mh2;
  if (hpara == 0) Mt2 = 0.18;
  else Mt2 = hpara[6];
  if (Hpara == 0) Mh2 = 0.67;
  else Mh2 = Hpara[6];
  double kt2 = 0.25 * Mt2 / (Mt2 + 0.25);
  double pt2 = 0.20 * Mh2 / (Mh2 + 0.20);
  double Pt2 = pt2 + z*z*kt2;
  double h1t[6], H1[6];
  const double Mneutron = 0.93957;
  const double Mpion = 0.13957;
  double pre = 0.5 * z * pow(kt2,2) * pt2 * pow(Pt,3) / (pow(Mneutron,2) * Mpion * pow(Pt2,3));
  lst->h1tn(hvar, h1t, hpara);
  lst->H1pip(Hvar, H1, Hpara);
  FF[0] = pre * x * (pow(2.0/3.0,2) * (h1t[0]*H1[0] + h1t[3]*H1[3]) + pow(1.0/3.0,2) * (h1t[1]*H1[1] + h1t[2]*H1[2] + h1t[4]*H1[4] + h1t[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi+
  lst->H1pim(Hvar, H1, Hpara);
  FF[1] = pre * x * (pow(2.0/3.0,2) * (h1t[0]*H1[0] + h1t[3]*H1[3]) + pow(1.0/3.0,2) * (h1t[1]*H1[1] + h1t[2]*H1[2] + h1t[4]*H1[4] + h1t[5]*H1[5])) * exp(-pow(Pt,2)/Pt2) / (M_PI * Pt2);//pi-
  return 0;
}

int Lstructure::FUTsin3HmSN(double * AZ, double * var, double * FF, double * hpara = 0, double * Hpara = 0){
  //AZ: Np, Nn, Polp, Poln
  double Np = AZ[0];
  double Nn = AZ[1];
  double Polp = AZ[2];
  double Poln = AZ[3];
  double FFp[2], FFn[2];
  lst->FUTsin3HmSp(var, FFp, hpara, Hpara);
  lst->FUTsin3HmSn(var, FFn, hpara, Hpara);
  FF[0] = Np * Polp * FFp[0] + Nn * Poln * FFn[0];
  FF[1] = Np * Polp * FFp[1] + Nn * Poln * FFn[1];
  return 0;
}

int Lstructure::AsinHmSp(double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double FFT[2], FFU[2];
  lst->FUTsinHmSp(var, FFT, Tf, TD);
  lst->FUUTp(var, FFU, Uf, UD);
  Asym[0] = FFT[0] / FFU[0];//pi+
  Asym[1] = FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::AsinHmSn(double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double FFT[2], FFU[2];
  lst->FUTsinHmSn(var, FFT, Tf, TD);
  lst->FUUTn(var, FFU, Uf, UD);
  Asym[0] = FFT[0] / FFU[0];//pi+
  Asym[1] = FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::AsinHmSN(double * AZ, double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double FFT[2], FFU[2];
  lst->FUTsinHmSN(AZ, var, FFT, Tf, TD);
  lst->FUUTN(AZ, var, FFU, Uf, UD);
  Asym[0] = FFT[0] / FFU[0];//pi+
  Asym[1] = FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::AsinHpSp(double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double y = kin[1];
  const double Mnucleon = 0.939;
  double gamma = 2.0 * Mnucleon * var[0] / sqrt(var[1]);
  double eps = (1.0 - y - 0.25*pow(y*gamma,2)) / (1.0 - y + 0.5*y*y + 0.25*pow(y*gamma,2));
  double FFT[2], FFU[2];
  lst->FUTsinHpSp(var, FFT, Tf, TD);
  lst->FUUTp(var, FFU, Uf, UD);
  Asym[0] = eps * FFT[0] / FFU[0];//pi+
  Asym[1] = eps * FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::AsinHpSn(double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double y = kin[1];
  const double Mnucleon = 0.939;
  double gamma = 2.0 * Mnucleon * var[0] / sqrt(var[1]);
  double eps = (1.0 - y - 0.25*pow(y*gamma,2)) / (1.0 - y + 0.5*y*y + 0.25*pow(y*gamma,2));
  double FFT[2], FFU[2];
  lst->FUTsinHpSn(var, FFT, Tf, TD);
  lst->FUUTn(var, FFU, Uf, UD);
  Asym[0] = eps * FFT[0] / FFU[0];//pi+
  Asym[1] = eps * FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::AsinHpSN(double * AZ, double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double y = kin[1];
  const double Mnucleon = 0.939;
  double gamma = 2.0 * Mnucleon * var[0] / sqrt(var[1]);
  double eps = (1.0 - y - 0.25*pow(y*gamma,2)) / (1.0 - y + 0.5*y*y + 0.25*pow(y*gamma,2));
  double FFT[2], FFU[2];
  lst->FUTsinHpSN(AZ, var, FFT, Tf, TD);
  lst->FUUTN(AZ, var, FFU, Uf, UD);
  Asym[0] = eps * FFT[0] / FFU[0];//pi+
  Asym[1] = eps * FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::Asin3HmSp(double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double y = kin[1];
  const double Mnucleon = 0.939;
  double gamma = 2.0 * Mnucleon * var[0] / sqrt(var[1]);
  double eps = (1.0 - y - 0.25*pow(y*gamma,2)) / (1.0 - y + 0.5*y*y + 0.25*pow(y*gamma,2));
  double FFT[2], FFU[2];
  lst->FUTsin3HmSp(var, FFT, Tf, TD);
  lst->FUUTp(var, FFU, Uf, UD);
  Asym[0] = eps * FFT[0] / FFU[0];//pi+
  Asym[1] = eps * FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::Asin3HmSn(double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double y = kin[1];
  const double Mnucleon = 0.939;
  double gamma = 2.0 * Mnucleon * var[0] / sqrt(var[1]);
  double eps = (1.0 - y - 0.25*pow(y*gamma,2)) / (1.0 - y + 0.5*y*y + 0.25*pow(y*gamma,2));
  double FFT[2], FFU[2];
  lst->FUTsin3HmSn(var, FFT, Tf, TD);
  lst->FUUTn(var, FFU, Uf, UD);
  Asym[0] = eps * FFT[0] / FFU[0];//pi+
  Asym[1] = eps * FFT[1] / FFU[1];//pi-
  return 0;
}

int Lstructure::Asin3HmSN(double * AZ, double * kin, double * Asym, double * Tf = 0, double * TD = 0, double * Uf = 0, double * UD = 0){
  //kin: x, y, z, Q2, Pt
  double var[4] = {kin[0], kin[3], kin[2], kin[4]};
  double y = kin[1];
  const double Mnucleon = 0.939;
  double gamma = 2.0 * Mnucleon * var[0] / sqrt(var[1]);
  double eps = (1.0 - y - 0.25*pow(y*gamma,2)) / (1.0 - y + 0.5*y*y + 0.25*pow(y*gamma,2));
  double FFT[2], FFU[2];
  lst->FUTsin3HmSN(AZ, var, FFT, Tf, TD);
  lst->FUUTN(AZ, var, FFU, Uf, UD);
  Asym[0] = eps * FFT[0] / FFU[0];//pi+
  Asym[1] = eps * FFT[1] / FFU[1];//pi-
  return 0;
}

double Lstructure::xsprefactor(double * kin){
  //kin: x, y, z, Q2, Pt, phih, phiS
  double x = kin[0];
  double y = kin[1];
  double Q2 = kin[3];
  const double alpha = 1.0/137.036;
  const double Mnucleon = 0.939;
  double gamma = 2.0 * x * Mnucleon / sqrt(Q2);
  double epsilon = (1.0 - y - 0.25 * pow(y*gamma,2)) / (1.0 - y + 0.5*y*y + 0.25*pow(y*gamma,2));
  return (alpha*alpha*y*y) / (2.0*(1.0-epsilon)*x*y*Q2) * (1.0 + gamma*gamma/2.0/x);//Ref. JHEP 02 (2007) 093
}

int Lstructure::CalcVariables(double * lab, double * phys){
  //lab: Ebeam, p_ele, theta_ele, phi_ele, p_pion, theta_pion, phi_pion
  const double degtorad = M_PI / 180.0;
  double Ebeam = lab[0];
  double p_ele = lab[1];
  double theta_ele = lab[2] * degtorad;
  double phi_ele = lab[3] * degtorad;
  double p_pion = lab[4];
  double theta_pion = lab[5] * degtorad;
  double phi_pion = lab[6] * degtorad;
  const double Mnucleon = 0.939;
  const double Mpion = 0.13957;
  double q0 = Ebeam - p_ele;
  double q1 = - p_ele * sin(theta_ele);
  double q3 = Ebeam - p_ele * cos(theta_ele);
  double qt = sqrt(q1*q1 + q3*q3);
  double Q2 = 2.0 * Ebeam * p_ele * (1.0 - cos(theta_ele));
  double x = Q2 / (2.0 * Mnucleon * q0);
  double y = q0 / Ebeam;
  double Epion = sqrt(Mpion * Mpion + p_pion * p_pion);
  double ppion1 = p_pion * sin(theta_pion) * cos(phi_pion - phi_ele);
  double ppion2 = p_pion * sin(theta_pion) * sin(phi_pion - phi_ele);
  double ppion3 = p_pion * cos(theta_pion);
  double z = Epion / q0;
  double W2 = Mnucleon * Mnucleon + (1.0 - x) / x * Q2;
  double Wp2 = W2 + Mpion * Mpion - 2.0 * ((Mnucleon + q0) * Epion - q1 * ppion1 - q3 * ppion3);
  double ct = q3 / qt;
  double st = q1 / qt;
  //double ppionz = ppion3 * ct + ppion1 * st;
  double ppionx = -ppion3 * st + ppion1 * ct;
  double Pt = sqrt(ppionx * ppionx + ppion2 * ppion2);
  double phi_h;
  if (ppion2 >= 0) phi_h = acos(ppionx / Pt);
  else phi_h = -acos(ppionx / Pt);
  double phi_S = 2.0 * M_PI - phi_ele;
  //phys: x, y, z, Q2, Pt, phih, phiS, W, Wp
  phys[0] = x;
  phys[1] = y;
  phys[2] = z;
  phys[3] = Q2;
  phys[4] = Pt;
  phys[5] = phi_h;
  phys[6] = phi_S;
  phys[7] = sqrt(W2);
  phys[8] = sqrt(Wp2);
  return 0;
}

double Lstructure::Jacobian(double * lab){
  //lab: Ebeam, p_ele, theta_ele, phi_ele, p_pion, theta_pion, phi_pion
  const double degtorad = M_PI / 180.0;
  const double Mnucleon = 0.939;
  const double Mpion = 0.13957;
  double Ebeam = lab[0];
  double p_ele = lab[1];
  double theta_ele = lab[2] * degtorad;
  double phi_ele = lab[3] * degtorad;
  double p_pion = lab[4];
  double theta_pion = lab[5] * degtorad;
  double phi_pion = lab[6] * degtorad;
  double p3 = Ebeam - p_ele * cos(theta_ele);
  double pT = p_pion * sqrt(1.0 + (-pow(p3 * cos(theta_pion), 2) + p3 * p_ele * sin(theta_ele) * sin(2.0 * theta_pion) * cos(phi_ele - phi_pion) - pow(p_ele * sin(theta_ele) * sin(theta_pion) * cos(phi_ele - phi_pion), 2)) / (pow(p3, 2) + pow(p_ele * sin(theta_ele), 2)) );
  double factor = (p_pion * p_ele * sin(theta_ele)) / (Mnucleon * pow(Ebeam - p_ele, 2) * sqrt(pow(p_pion, 2) + pow(Mpion, 2)));
  double jnum = 4.0 * p_pion * sin(theta_pion) * (p_ele * sin(theta_ele) * sin(theta_pion) * cos(phi_ele - phi_pion) - p3 * cos(theta_pion));
  double jden = sqrt(std::abs( -8.0 * pow(p3, 2) * cos(2.0 * theta_pion) + p_ele * (16.0 * p3 * sin(theta_ele) * sin(2.0 * theta_pion) * cos(phi_ele - phi_pion) - p_ele * (8.0 * pow(sin(theta_ele) * sin(theta_pion), 2) * cos(2.0 * (phi_ele - phi_pion)) + cos(2.0 * (theta_ele - theta_pion)) + cos(2.0 * (theta_ele + theta_pion)) + 6.0 * cos(2.0 * theta_ele) - 2.0 * cos(2.0 * theta_pion))) + 8.0 * pow(p3, 2) + 6.0 * pow(p_ele, 2)));
  return std::abs((2.0 * pT * factor * jnum) / (jden * sin(theta_ele) * sin(theta_pion)));
}

int Lstructure::sigmaUUT(double * AZ, double * lab, double * xs, double * fpara = 0, double * Dpara = 0){
  double phys[9];// x, y, z, Q2, Pt, phih, phiS, W, Wp
  lst->CalcVariables(lab, phys);
  double var[4] = {phys[0], phys[3], phys[2], phys[4]};
  double FF[2];
  lst->FUUTN(AZ, var, FF, fpara, Dpara);
  double pre = lst->xsprefactor(phys);
  double jac = lst->Jacobian(lab);
  xs[0] = jac * pre * FF[0];//pi+
  xs[1] = jac * pre * FF[1];//pi-
  return 0;
}

#endif
