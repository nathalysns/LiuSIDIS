#include "splitsusy.h"

int main(const int argc, const char * argv[]){

  
  cout << d1u(200.0, 200.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;
  cout << d1d(200.0, 200.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;

  cout << d2u(200.0, 200.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;
  cout << d2d(200.0, 200.0, 1.025 * sqrt(2.0) * MW / VEV, 1.025 * sqrt(2.0) * MW / VEV) << endl;

  return 0;
}
