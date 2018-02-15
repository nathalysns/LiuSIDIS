#ifndef _F2_H_
#define _F2_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "LHAPDF/LHAPDF.h"

const LHAPDF::PDF * xpdf;

double F2(const double x, const double Q2){
  double result = pow(2.0 / 3.0, 2) * (xpdf->xfxQ2(2, x, Q2) + xpdf->xfxQ2(-2, x, Q2) + xpdf->xfxQ2(4, x, Q2) + xpdf->xfxQ2(-4, x, Q2))
    + pow(1.0 / 3.0, 2) * (xpdf->xfxQ2(1, x, Q2) + xpdf->xfxQ2(-1, x, Q2) + xpdf->xfxQ2(3, x, Q2) + xpdf->xfxQ2(-3, x, Q2) + xpdf->xfxQ2(5, x, Q2) + xpdf->xfxQ2(-5, x, Q2));
  return result;
}


#endif
