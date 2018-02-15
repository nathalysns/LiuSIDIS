#include "F2.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2) {
    return 0;
  }

  const LHAPDF::PDF * xf[101];
  for (int i = 0; i <= 100; i++){
    xf[i] = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", i);
  }

  ifstream infile(argv[1]);

  double x, Q2;
  double central, sum, error;

  FILE * fs = fopen("out_NNPDF.dat", "w");
  fprintf(fs, "x\tQ2\tvalue\terror\n");

  while (infile >> x >> Q2){
    xpdf = xf[0];
    central = F2(x, Q2);
    sum = 0;
    for (int i = 1; i <= 100; i++){
      xpdf = xf[i];
      sum += pow(F2(x, Q2) - central, 2);
    }
    error = sqrt(sum / 100.0);
    fprintf(fs, "%.3E\t%.3E\t%.3E\t%.3E\n",
	    x, Q2, central, error);
  }
  infile.close();

  fclose(fs);

  return 0;
}
