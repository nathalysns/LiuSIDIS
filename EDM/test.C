#include "splitsusy.h"

int main(const int argc, const char * argv[]){


  cout << d1d(600.0, 190.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;
  cout << d2d(600.0, 190.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;
  cout << d3d(600.0, 190.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;

  cout << d1d(620.0, 190.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;
  cout << d2d(620.0, 190.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;
  cout << d3d(620.0, 190.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;

  return 0;
}
