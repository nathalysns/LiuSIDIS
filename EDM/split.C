#include <iostream>
#include <fstream>
#include <cmath>

#include "gsl/gsl_sf_dilog.h"

using namespace std;


const double MH = 125.09;//higgs mass
const double MZ = 91.1876;//Z mass
const double MW = 80.385;//W mass
const double sW2 = 0.23155;//mixing angle

double F_gH(const double rho, const double R);
double F_ZH(const double r, const double rho, const double R);
double F_WW1(const double rho, const double R);
double F_WW2(const double rho, const double R);
double A_ZH(const double rho);
double B_ZH(const double r, const double rho);
double A_WW1(const double rho);
double B_WW1(const double rho);
double A_WW2(const double rho);
double B_WW2(const double rho);

double Li2(const double z){//real part of the dilog function
  return gsl_sf_dilog(z);
}

int main(const int argc, const char * argv[]){

  cout << Li2(1.5) << endl;

  return 0;
}



/////
double F_gH(const double rho, const double R){
  return -0.5 * log(R) - 1.0 + (rho + 1.0) * log(rho) / (4.0 * (rho - 1.0));
}

double F_ZH(const double r, const double rho, const double R){
  return A_ZH(rho) * log(R) + B_ZH(r, rho);
}

double F_WW1(const double rho, const double R){
  return A_WW1(rho) * log(R) + B_WW1(rho);
}

double F_WW2(const double rho, const double R){
  return A_WW2(rho) * log(R) + B_WW2(rho);
}

double A_ZH(const double rho){
  return ((rho - 1.0) * (2.0 - rho) - rho * log(rho)) / (4.0 * pow(rho - 1.0, 2)) + sW2 / 2.0;
}

double B_ZH(const double r, const double rho){
  return 0;
}

double A_WW1(const double rho){
  return 0;
}

double B_WW1(const double rho){
  return 0;
}

double A_WW2(const double rho){
  return 0;
}

double B_WW2(const double rho){
  return 0;
}
