#ifndef _SSUSY_H_
#define _SSUSY_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"

#include "gsl/gsl_sf_dilog.h"

double jloop(const double r){
  if (r < 0) return 0;
  if (r == 1.0) return 1.0;
  return r * log(r) / (r - 1.0);
}

double jloop(const double r, const double s){
  if (r < 0 || s < 0) return 0;
  if (r == s){
    if (s == 1.0) return 0.5;
    else return (s - 1.0 - log(s)) / pow(s - 1.0, 2);
  }
  return (jloop(r) - jloop(s)) / (r - s);
}







#endif
